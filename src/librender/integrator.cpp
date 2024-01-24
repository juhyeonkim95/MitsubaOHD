/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/renderproc.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>


MTS_NAMESPACE_BEGIN

Integrator::Integrator(const Properties &props)
 : NetworkedObject(props) { }

Integrator::Integrator(Stream *stream, InstanceManager *manager)
 : NetworkedObject(stream, manager) { }

bool Integrator::preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) { return true; }
void Integrator::postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) { }
void Integrator::serialize(Stream *stream, InstanceManager *manager) const {
    NetworkedObject::serialize(stream, manager);
}
void Integrator::configureSampler(const Scene *scene, Sampler *sampler) {
    /* Prepare the sampler for bucket-based rendering */
    sampler->setFilmResolution(scene->getFilm()->getCropSize(),
        getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)));
}
const Integrator *Integrator::getSubIntegrator(int idx) const { return NULL; }

SamplingIntegrator::SamplingIntegrator(const Properties &props)
 : Integrator(props) {

}

void SamplingIntegrator::addChild(const std::string &name, ConfigurableObject *child) {

}

SamplingIntegrator::SamplingIntegrator(Stream *stream, InstanceManager *manager)
 : Integrator(stream, manager) { }

void SamplingIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
    Integrator::serialize(stream, manager);
}

Spectrum SamplingIntegrator::E(const Scene *scene, const Intersection &its,
        const Medium *medium, Sampler *sampler, int nSamples, bool handleIndirect) const {
    Spectrum E(0.0f);
    RadianceQueryRecord query(scene, sampler);
    DirectSamplingRecord dRec(its);
    Frame frame(its.shFrame.n);

    sampler->generate(Point2i(0));
    for (int i=0; i<nSamples; i++) {
        /* Sample the direct illumination component */
        int maxIntermediateInteractions = -1;
        Spectrum directRadiance = scene->sampleAttenuatedEmitterDirect(
            dRec, its, medium, maxIntermediateInteractions, query.nextSample2D());

        if (!directRadiance.isZero()) {
            Float dp = dot(dRec.d, its.shFrame.n);
            if (dp > 0)
                E += directRadiance * dp;
        }

        /* Sample the indirect illumination component */
        if (handleIndirect) {
            query.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);
            Vector d = frame.toWorld(warp::squareToCosineHemisphere(query.nextSample2D()));
            ++query.depth;
            query.medium = medium;
            E += Li(RayDifferential(its.p, d, its.time), query) * M_PI;
        }

        sampler->advance();
    }

    return E / (Float) nSamples;
}

void SamplingIntegrator::cancel() {
    if (m_process)
        Scheduler::getInstance()->cancel(m_process);
}

bool SamplingIntegrator::render(Scene *scene,
        RenderQueue *queue, const RenderJob *job,
        int sceneResID, int sensorResID, int samplerResID) {
    ref<Scheduler> sched = Scheduler::getInstance();
    ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
    ref<Film> film = sensor->getFilm();

    size_t nCores = sched->getCoreCount();
    const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
    size_t sampleCount = sampler->getSampleCount();

    Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
        " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
        sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
        nCores == 1 ? "core" : "cores");

    /* This is a sampling-based integrator - parallelize */
    // ref<ParallelProcess>
    ref<BlockedRenderProcess> proc = new BlockedRenderProcess(job,
        queue, scene->getBlockSize());

    proc->setPixelFormat(
        film->getFrames() > 1 ? Bitmap::EMultiSpectrumAlphaWeight : Bitmap::ESpectrumAlphaWeight,
        (int) (film->getFrames() * SPECTRUM_SAMPLES + 2), false);
        
    int integratorResID = sched->registerResource(this);
    proc->bindResource("integrator", integratorResID);
    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);
    scene->bindUsedResources(proc);
    bindUsedResources(proc);
    sched->schedule(proc);

    m_process = proc;
    sched->wait(proc);
    m_process = NULL;
    sched->unregisterResource(integratorResID);

    return proc->getReturnStatus() == ParallelProcess::ESuccess;
}

void SamplingIntegrator::bindUsedResources(ParallelProcess *) const {
    /* Do nothing by default */
}

void SamplingIntegrator::wakeup(ConfigurableObject *parent,
    std::map<std::string, SerializableObject *> &) {
    /* Do nothing by default */
}

void SamplingIntegrator::renderBlock(const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {
    // printf("\nrenderBlock: %d, %d, %d\n", block->getBitmap()->getWidth(), block->getBitmap()->getHeight(), block->getBitmap()->getChannelCount());
        
    Float diffScaleFactor = 1.0f /
        std::sqrt((Float) sampler->getSampleCount());

    bool needsApertureSample = sensor->needsApertureSample();
    bool needsTimeSample = sensor->needsTimeSample();

    RadianceQueryRecord rRec(scene, sampler);
    Point2 apertureSample(0.5f);
    Float timeSample = 0.5f;
    RayDifferential sensorRay;

    block->clear();

    uint32_t queryType = RadianceQueryRecord::ESensorRay;

    if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
        queryType &= ~RadianceQueryRecord::EOpacity;

    for (size_t i = 0; i<points.size(); ++i) {
        Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
        if (stop)
            break;

        sampler->generate(offset);

        for (size_t j = 0; j<sampler->getSampleCount(); j++) {
            rRec.newQuery(queryType, sensor->getMedium());
            Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

            if (needsApertureSample)
                apertureSample = rRec.nextSample2D();
            if (needsTimeSample)
                timeSample = rRec.nextSample1D();

            Spectrum spec = sensor->sampleRayDifferential(
                sensorRay, samplePos, apertureSample, timeSample);

            sensorRay.scaleDifferential(diffScaleFactor);

            spec *= Li(sensorRay, rRec);
            block->put(samplePos, spec, rRec.alpha);
            sampler->advance();
        }
    }
}


MonteCarloIntegrator::MonteCarloIntegrator(const Properties &props) : SamplingIntegrator(props) {
    /* Depth to begin using russian roulette */
    m_rrDepth = props.getInteger("rrDepth", 5);

    /* Longest visualized path depth (\c -1 = infinite).
       A value of \c 1 will visualize only directly visible light sources.
       \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
    m_maxDepth = props.getInteger("maxDepth", -1);

    /**
     * This parameter specifies the action to be taken when the geometric
     * and shading normals of a surface don't agree on whether a ray is on
     * the front or back-side of a surface.
     *
     * When \c strictNormals is set to \c false, the shading normal has
     * precedence, and rendering proceeds normally at the risk of
     * introducing small light leaks (this is the default).
     *
     * When \c strictNormals is set to \c true, the random walk is
     * terminated when encountering such a situation. This may
     * lead to silhouette darkening on badly tesselated meshes.
     */
    m_strictNormals = props.getBoolean("strictNormals", false);

    /**
     * When this flag is set to true, contributions from directly
     * visible emitters will not be included in the rendered image
     */
    m_hideEmitters = props.getBoolean("hideEmitters", false);

    if (m_rrDepth <= 0)
        Log(EError, "'rrDepth' must be set to a value greater than zero!");

    if (m_maxDepth <= 0 && m_maxDepth != -1)
        Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
}

MonteCarloIntegrator::MonteCarloIntegrator(Stream *stream, InstanceManager *manager)
    : SamplingIntegrator(stream, manager) {
    m_rrDepth = stream->readInt();
    m_maxDepth = stream->readInt();
    m_strictNormals = stream->readBool();
    m_hideEmitters = stream->readBool();
}

void MonteCarloIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
    SamplingIntegrator::serialize(stream, manager);
    stream->writeInt(m_rrDepth);
    stream->writeInt(m_maxDepth);
    stream->writeBool(m_strictNormals);
    stream->writeBool(m_hideEmitters);
}

// std::vector<std::tuple<Spectrum, Float>> ExhaustivePathIntegrator::Li_helper(const RayDifferential &r, RadianceQueryRecord &rRec) const {
//     /* Some aliases and local variables */
//     const Scene *scene = rRec.scene;
//     Intersection &its = rRec.its;
//     RayDifferential ray(r);
//     Spectrum Li(0.0f);
//     bool scattered = false;

//     std::vector<std::tuple<Spectrum, Float>> results;

//     /* Perform the first ray intersection (or ignore if the
//         intersection has already been provided). */
//     rRec.rayIntersect(ray);
//     ray.mint = Epsilon;

//     Spectrum throughput(1.0f);
//     Float eta = 1.0f;
//     Float path_length = 0.0f;
//     path_length += its.t * eta;

//     Intersection first_its = rRec.its;
//     if (!first_its.isValid()) {
//         return results;
//     }
//     const BSDF *first_bsdf = first_its.getBSDF(ray);
//     Float first_path_length = path_length;

//     while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
//         if (!its.isValid()) {
//             /* If no intersection could be found, potentially return
//                 radiance from a environment luminaire if it exists */
//             // if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
//             //     && (!m_hideEmitters || scattered))
//             //     Spectrum Li = throughput * scene->evalEnvironment(ray);
//             break;
//         }

//         const BSDF *bsdf = its.getBSDF(ray);

//         /* Possibly include emitted radiance if requested */
//         if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
//             && (!m_hideEmitters || scattered))
//         {
//             Spectrum Li = throughput * its.Le(-ray.d);
//             results.push_back({Li, path_length});
//         }
//             // Li += throughput * its.Le(-ray.d);

//         /* Include radiance from a subsurface scattering model if requested */
//         //if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
//         //    Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

//         if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
//             || (m_strictNormals && dot(ray.d, its.geoFrame.n)
//                 * Frame::cosTheta(its.wi) >= 0)) {

//             /* Only continue if:
//                 1. The current path length is below the specifed maximum
//                 2. If 'strictNormals'=true, when the geometric and shading
//                     normals classify the incident direction to the same side */
//             break;
//         }

//         /* ==================================================================== */
//         /*                     Direct illumination sampling                     */
//         /* ==================================================================== */

//         /* Estimate the direct illumination if this is requested */
//         DirectSamplingRecord dRec(its);

//         if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
//             (bsdf->getType() & BSDF::ESmooth)) {

//             if(m_use_collimated_light && rRec.depth > 1){
//                 Point origin = first_its.p;
//                 Vector direction = (origin - its.p);
//                 Float distance = direction.length();
//                 Float invDist = 1.0f / distance;
//                 direction *= invDist;
//                 Ray ray_temp(origin, direction, Epsilon,
//                     distance*(1-ShadowEpsilon), ray.time);
//                 Intersection its_temp;

//                 if(!scene->rayIntersect(ray_temp, its_temp)){
//                     Spectrum value(1.0);

//                     /* Allocate a record for querying the BSDF */
//                     BSDFSamplingRecord bRec(its, its.toLocal(direction), ERadiance);

//                     /* Evaluate BSDF * cos(theta) */
//                     const Spectrum bsdfVal = bsdf->eval(bRec);

//                     /* Prevent light leaks due to the use of shading normals */
//                     if (!bsdfVal.isZero() && (!m_strictNormals
//                             || dot(its.geoFrame.n, direction) * Frame::cosTheta(bRec.wo) > 0)) {

//                         /* Calculate prob. of having generated that direction
//                         using BSDF sampling */
//                         Float bsdfPdf = 0;

//                         /* Weight using the power heuristic */
//                         Float weight = 1.0; //miWeight(dRec.pdf, bsdfPdf);

//                         Float em_path_length = path_length + distance + first_path_length;
                    
//                         /* Allocate a record for querying the BSDF */
//                         BSDFSamplingRecord bRec2(first_its, first_its.toLocal(-direction), ERadiance);

//                         /* Evaluate BSDF * cos(theta) */
//                         const Spectrum bsdfVal2 = first_bsdf->eval(bRec2);
                        
//                         Spectrum Li = throughput * value * bsdfVal * weight * bsdfVal2;
                        
//                         results.push_back({Li, em_path_length});

//                         // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
//                     }
//                 }
//             } else if (m_use_collimated_light && rRec.depth == 1) {
//                 Point origin = ray.o;
//                 Vector direction = (origin - its.p);
//                 Float distance = direction.length();
//                 Float invDist = 1.0f / distance;
//                 direction *= invDist;

//                 Spectrum value(1.0);

//                 /* Allocate a record for querying the BSDF */
//                 BSDFSamplingRecord bRec(its, its.toLocal(direction), ERadiance);

//                 /* Evaluate BSDF * cos(theta) */
//                 const Spectrum bsdfVal = bsdf->eval(bRec);

//                 /* Prevent light leaks due to the use of shading normals */
//                 if (!bsdfVal.isZero() && (!m_strictNormals
//                         || dot(its.geoFrame.n, direction) * Frame::cosTheta(bRec.wo) > 0)) {

//                     /* Calculate prob. of having generated that direction
//                     using BSDF sampling */
//                     Float bsdfPdf = 0;

//                     /* Weight using the power heuristic */
//                     Float weight = 1.0; //miWeight(dRec.pdf, bsdfPdf);

//                     Float em_path_length = path_length + distance;
                    
//                     Spectrum Li = throughput * value * bsdfVal * weight;
//                     results.push_back({Li, em_path_length});
//                 }
                
//             } else {

//             Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
//             if (!value.isZero()) {
//                 const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

//                 /* Allocate a record for querying the BSDF */
//                 BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

//                 /* Evaluate BSDF * cos(theta) */
//                 const Spectrum bsdfVal = bsdf->eval(bRec);

//                 /* Prevent light leaks due to the use of shading normals */
//                 if (!bsdfVal.isZero() && (!m_strictNormals
//                         || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

//                     /* Calculate prob. of having generated that direction
//                         using BSDF sampling */
//                     Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
//                         ? bsdf->pdf(bRec) : 0;

//                     /* Weight using the power heuristic */
//                     Float weight = miWeight(dRec.pdf, bsdfPdf);

//                     Float em_path_length = path_length + dRec.dist;
//                     Spectrum Li = throughput * value * bsdfVal * weight;
                    
//                     results.push_back({Li, em_path_length});

//                     // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
//                 }
//             }}
//         }

//         /* ==================================================================== */
//         /*                            BSDF sampling                             */
//         /* ==================================================================== */

//         /* Sample BSDF * cos(theta) */
//         Float bsdfPdf;
//         BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
//         Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
//         if (bsdfWeight.isZero())
//             break;

//         scattered |= bRec.sampledType != BSDF::ENull;

//         /* Prevent light leaks due to the use of shading normals */
//         const Vector wo = its.toWorld(bRec.wo);
//         Float woDotGeoN = dot(its.geoFrame.n, wo);
//         if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
//             break;

//         bool hitEmitter = false;
//         Spectrum value;

//         /* Trace a ray in this direction */
//         ray = Ray(its.p, wo, ray.time);
//         if (scene->rayIntersect(ray, its)) {
//             /* Intersected something - check if it was a luminaire */
//             if (its.isEmitter()) {
//                 value = its.Le(-ray.d);
//                 dRec.setQuery(ray, its);
//                 hitEmitter = true;
//             }
//         } else {
//             /* Intersected nothing -- perhaps there is an environment map? */
//             const Emitter *env = scene->getEnvironmentEmitter();

//             if (env) {
//                 if (m_hideEmitters && !scattered)
//                     break;

//                 value = env->evalEnvironment(ray);
//                 if (!env->fillDirectSamplingRecord(dRec, ray))
//                     break;
//                 hitEmitter = true;
//             } else {
//                 break;
//             }
//         }

//         /* Keep track of the throughput and relative
//             refractive index along the path */
//         throughput *= bsdfWeight;
//         eta *= bRec.eta;
//         path_length += its.t * eta;

//         /* If a luminaire was hit, estimate the local illumination and
//             weight using the power heuristic */
//         if (hitEmitter &&
//             (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
//             /* Compute the prob. of generating that direction using the
//                 implemented direct illumination sampling technique */
//             const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
//                 scene->pdfEmitterDirect(dRec) : 0;
            
//             Spectrum Li = throughput * value * miWeight(bsdfPdf, lumPdf);
//             results.push_back({Li, path_length});
//             // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
//         }

//         /* ==================================================================== */
//         /*                         Indirect illumination                        */
//         /* ==================================================================== */

//         /* Set the recursive query type. Stop if no surface was hit by the
//             BSDF sample or if indirect illumination was not requested */
//         if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
//             break;
//         rRec.type = RadianceQueryRecord::ERadianceNoEmission;

//         if (rRec.depth++ >= m_rrDepth) {
//             /* Russian roulette: try to keep path weights equal to one,
//                 while accounting for the solid angle compression at refractive
//                 index boundaries. Stop with at least some probability to avoid
//                 getting stuck (e.g. due to total internal reflection) */

//             Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
//             if (rRec.nextSample1D() >= q)
//                 break;
//             throughput /= q;
//         }
//     }

//     /* Store statistics */
//     avgPathLength.incrementBase();
//     avgPathLength += rRec.depth;

//     return results;
// }

std::string RadianceQueryRecord::toString() const {
    std::ostringstream oss;
    oss << "RadianceQueryRecord[" << endl
        << "  type = { ";
    if (type & EEmittedRadiance) oss << "emitted ";
    if (type & ESubsurfaceRadiance) oss << "subsurface ";
    if (type & EDirectSurfaceRadiance) oss << "direct ";
    if (type & EIndirectSurfaceRadiance) oss << "indirect ";
    if (type & ECausticRadiance) oss << "caustic ";
    if (type & EDirectMediumRadiance) oss << "inscatteredDirect ";
    if (type & EIndirectMediumRadiance) oss << "inscatteredIndirect ";
    if (type & EDistance) oss << "distance ";
    if (type & EOpacity) oss << "opacity ";
    if (type & EIntersection) oss << "intersection ";
    oss << "}," << endl
        << "  depth = " << depth << "," << endl
        << "  its = " << indent(its.toString()) << endl
        << "  alpha = " << alpha << "," << endl
        << "  extra = " << extra << "," << endl
        << "]" << endl;
    return oss.str();
}


MTS_IMPLEMENT_CLASS(Integrator, true, NetworkedObject)
MTS_IMPLEMENT_CLASS(SamplingIntegrator, true, Integrator)
MTS_IMPLEMENT_CLASS(MonteCarloIntegrator, true, SamplingIntegrator)
// MTS_IMPLEMENT_CLASS(ExhaustivePathIntegrator, true, MonteCarloIntegrator)
MTS_NAMESPACE_END
