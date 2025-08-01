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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include "fmcw_interface.h"

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

/*! \plugin{path}{Path tracer}
 * \order{2}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *        which the implementation will start to use the ``russian roulette''
 *        path termination criterion. \default{\code{5}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See the description below
 *        for details.\default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 *
 * This integrator implements a basic path tracer and is a \emph{good default choice}
 * when there is no strong reason to prefer another method.
 *
 * To use the path tracer appropriately, it is instructive to know roughly how
 * it works: its main operation is to trace many light paths using \emph{random walks}
 * starting from the sensor. A single random walk is shown below, which entails
 * casting a ray associated with a pixel in the output image and searching for
 * the first visible intersection. A new direction is then chosen at the intersection,
 * and the ray-casting step repeats over and over again (until one of several
 * stopping criteria applies).
 * \begin{center}
 * \includegraphics[width=.7\textwidth]{images/integrator_path_figure.pdf}
 * \end{center}
 * At every intersection, the path tracer tries to create a connection to
 * the light source in an attempt to find a \emph{complete} path along which
 * light can flow from the emitter to the sensor. This of course only works
 * when there is no occluding object between the intersection and the emitter.
 *
 * This directly translates into a category of scenes where
 * a path tracer can be expected to produce reasonable results: this is the case
 * when the emitters are easily ``accessible'' by the contents of the scene. For instance,
 * an interior scene that is lit by an area light will be considerably harder
 * to render when this area light is inside a glass enclosure (which
 * effectively counts as an occluder).
 *
 * Like the \pluginref{direct} plugin, the path tracer internally relies on multiple importance
 * sampling to combine BSDF and emitter samples. The main difference in comparison
 * to the former plugin is that it considers light paths of arbitrary length to compute
 * both direct and indirect illumination.
 *
 * For good results, combine the path tracer with one of the
 * low-discrepancy sample generators (i.e. \pluginref{ldsampler},
 * \pluginref{halton}, or \pluginref{sobol}).
 *
 * \paragraph{Strict normals:}\label{sec:strictnormals}
 * Triangle meshes often rely on interpolated shading normals
 * to suppress the inherently faceted appearance of the underlying geometry. These
 * ``fake'' normals are not without problems, however. They can lead to paradoxical
 * situations where a light ray impinges on an object from a direction that is classified as ``outside''
 * according to the shading normal, and ``inside'' according to the true geometric normal.
 *
 * The \code{strictNormals}
 * parameter specifies the intended behavior when such cases arise. The default (\code{false}, i.e. ``carry on'')
 * gives precedence to information given by the shading normal and considers such light paths to be valid.
 * This can theoretically cause light ``leaks'' through boundaries, but it is not much of a problem in practice.
 *
 * When set to \code{true}, the path tracer detects inconsistencies and ignores these paths. When objects
 * are poorly tesselated, this latter option may cause them to lose a significant amount of the incident
 * radiation (or, in other words, they will look dark).
 *
 * The bidirectional integrators in Mitsuba (\pluginref{bdpt}, \pluginref{pssmlt}, \pluginref{mlt} ...)
 * implicitly have \code{strictNormals} set to \code{true}. Hence, another use of this parameter
 * is to match renderings created by these methods.
 *
 * \remarks{
 *    \item This integrator does not handle participating media
 *    \item This integrator has poor convergence properties when rendering
 *    caustics and similar effects. In this case, \pluginref{bdpt} or
 *    one of the photon mappers may be preferable.
 * }
 */
class MIFMCWPSDPathTracer : public MonteCarloIntegrator, public FMCWInterface {
public:
    MIFMCWPSDPathTracer(const Properties &props)
        : MonteCarloIntegrator(props), FMCWInterface(props) {
        m_max_distance = props.getFloat("max_distance", 10);
        m_min_distance = props.getFloat("min_distance", 0);
        printf("[Transient INFO]\n");
        printf("\tmax_distance: %f\n", m_max_distance);
        printf("\tmin_distance: %f\n", m_min_distance);
    }

    /// Unserialize from a binary data stream
    MIFMCWPSDPathTracer(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        NotImplementedError("Li");
    }

    std::tuple<size_t, Float> freq_to_index(Float frequency) const {
        // frequency * 1e6 * T * 1e-6 * 3e8 / (B * 1e9)
        Float path_length = frequency * m_T / m_B * 3.0 / 10.0;
        
        if(path_length < m_min_distance || path_length > m_max_distance){
            printf("OUT OF RANGE!, %.2f\n", path_length);
            return {-1.0, -1.0};
        }
        Float idx = (path_length - m_min_distance) / (m_max_distance - m_min_distance) * m_M;
        size_t idx_int = (size_t)(idx);
        idx_int = std::clamp(idx_int, (size_t)0, m_M-1);
        Float float_idx = idx - idx_int;
        float_idx = std::clamp(float_idx, 0.0, 1.0);
        return {idx_int, float_idx};
    }

    std::vector<std::tuple<Spectrum, Float, Float, Float>> Li_helper(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        std::vector<std::tuple<Spectrum, Float, Float, Float>> results;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;
        Float path_length = 0.0f;
        Float path_velocity = 0.0f;
        Vector prev_velocity = scene->getSensor()->getVelocity();//Vector(0.0f);

        path_length += its.t * eta;
        // prev_velocity = its.getVelocity();

        Intersection first_its = rRec.its;
        if (!first_its.isValid()) {
            return results;
        }
        const BSDF *first_bsdf = first_its.getBSDF(ray);
        Float first_path_length = path_length;
        Vector first_velocity = its.getVelocity();
        Float first_path_velocity = dot( (first_velocity - prev_velocity), ray.d);
        Float path_pdf = 1.0;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                // if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                //     && (!m_hideEmitters || scattered))
                //     Spectrum Li = throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            Vector velocity = its.getVelocity();
            path_velocity += dot( (velocity - prev_velocity), ray.d);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
            {
                Spectrum Li = throughput * its.Le(-ray.d);
                results.push_back({Li, path_pdf, path_length, path_velocity});
            }
                // Li += throughput * its.Le(-ray.d);

            /* Include radiance from a subsurface scattering model if requested */
            //if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
            //    Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

            if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                    * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);

            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) {

                if(m_use_collimated_light && rRec.depth > 1){
                    Point origin = first_its.p;
                    Vector direction = (origin - its.p);
                    Float distance = direction.length();
                    Float invDist = 1.0f / distance;
                    direction *= invDist;
                    Ray ray_temp(origin, direction, Epsilon,
                        distance*(1-ShadowEpsilon), ray.time);
                    Intersection its_temp;

                    if(!scene->rayIntersect(ray_temp, its_temp)){
                        Spectrum value(1.0);

                        /* Allocate a record for querying the BSDF */
                        BSDFSamplingRecord bRec(its, its.toLocal(direction), ERadiance);

                        /* Evaluate BSDF * cos(theta) */
                        const Spectrum bsdfVal = bsdf->eval(bRec);

                        /* Prevent light leaks due to the use of shading normals */
                        if (!bsdfVal.isZero() && (!m_strictNormals
                                || dot(its.geoFrame.n, direction) * Frame::cosTheta(bRec.wo) > 0)) {

                            /* Calculate prob. of having generated that direction
                            using BSDF sampling */
                            Float bsdfPdf = 0;

                            /* Weight using the power heuristic */
                            Float weight = 1.0; //miWeight(dRec.pdf, bsdfPdf);

                            Float em_path_length = path_length + distance + first_path_length;
                            Float em_velocity = path_velocity + dot(first_velocity, direction) + first_path_velocity;
                        
                            /* Allocate a record for querying the BSDF */
                            BSDFSamplingRecord bRec2(first_its, first_its.toLocal(-direction), ERadiance);

                            /* Evaluate BSDF * cos(theta) */
                            const Spectrum bsdfVal2 = first_bsdf->eval(bRec2);
                            
                            Spectrum Li = throughput * value * bsdfVal * weight * bsdfVal2;
                            
                            results.push_back({Li, path_pdf, em_path_length, em_velocity});

                            // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
                        }
                    }
                } else if (m_use_collimated_light && rRec.depth == 1) {
                    Point origin = ray.o;
                    Vector direction = (origin - its.p);
                    Float distance = direction.length();
                    Float invDist = 1.0f / distance;
                    direction *= invDist;

                    Spectrum value(1.0);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(direction), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_strictNormals
                            || dot(its.geoFrame.n, direction) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                        using BSDF sampling */
                        Float bsdfPdf = 0;

                        /* Weight using the power heuristic */
                        Float weight = 1.0; //miWeight(dRec.pdf, bsdfPdf);

                        Float em_path_length = path_length + distance;
                        Float em_velocity = path_velocity + first_path_velocity;
                        
                        Spectrum Li = throughput * value * bsdfVal * weight;
                        results.push_back({Li, path_pdf, em_path_length, em_velocity});
                    }
                    
                } else {

                Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
                if (!value.isZero()) {
                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_strictNormals
                            || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                           using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                            ? bsdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        Float weight = miWeight(dRec.pdf, bsdfPdf);

                        Float em_path_length = path_length + dRec.dist;
                        Float em_velocity = path_velocity + dot(velocity, -dRec.d);
                        Spectrum Li = throughput * value * bsdfVal * weight;
                        
                        results.push_back({Li, path_pdf * dRec.pdf, em_path_length, em_velocity});

                        // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
                    }
                }}
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */
            throughput *= bsdfWeight;
            path_pdf *= bsdfPdf;
            eta *= bRec.eta;
            path_length += its.t * eta;
            prev_velocity = velocity * eta;

            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                    scene->pdfEmitterDirect(dRec) : 0;
                
                Spectrum Li = throughput * value * miWeight(bsdfPdf, lumPdf);
                results.push_back({Li, path_pdf, path_length, path_velocity});
                // put_Li_to_block(Li, em_path_length, block, sampled_pos, rRec.alpha);
            }

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
        }

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;

        return results;
    }


    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MIFMCWPSDPathTracer[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    // Overloaded
    void renderBlock(const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {
       
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

        Float *temp = (Float *) alloca(sizeof(Float) * (2 * m_M * SPECTRUM_SAMPLES + 2));

        for (size_t i = 0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            if (stop)
                break;

            sampler->generate(offset);


            // for each samples
            for (size_t j = 0; j<sampler->getSampleCount(); j++) {
                rRec.newQuery(queryType, sensor->getMedium());
                
                // get sample
                // if(use_collimated){
                //     m_laser_mrad 
                // }

                // Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

                Vector2 jitter = Vector2(rRec.nextSample2D()) * 2.0 - Vector2(1.0);
                jitter = jitter * m_fov_error + Vector2(0.5);
                Point2 samplePos = Point2(offset) + jitter;

                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                if (needsTimeSample)
                    timeSample = rRec.nextSample1D();

                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos, apertureSample, timeSample);

                sensorRay.scaleDifferential(diffScaleFactor);

                // all time samples
                std::vector<std::tuple<Spectrum, Float, Float, Float>> results = Li_helper(sensorRay, rRec);
                
                // put all time samples

                // (1) initialize array
                for(size_t k=0; k<(2 * m_M * SPECTRUM_SAMPLES + 2); k++){
                    temp[k] = 0.0;
                }

                // (2) for each sampled path -> put to result
                for (size_t k=0; k<results.size(); k++){
                    auto[Li, path_pdf, path_length, path_velocity] = results.at(k);
                    Li = Li * spec;

                    Float f_up = get_fmcw_weight_velocity_freq(path_length, path_velocity);
                    Float f_down = get_fmcw_weight_velocity_inverted_freq(path_length, path_velocity);
                    
                    // up & down chirp
                    auto[idx1, remainder1] = freq_to_index(f_up);
                    auto[idx2, remainder2] = freq_to_index(f_down);

                    size_t idx11 = std::clamp(idx1 + 1, (size_t)0, m_M-1);
                    size_t idx22 = std::clamp(idx2 + 1, (size_t)0, m_M-1);

                    // (3) put to each time samples
                    if(idx1 >= 0){
                        for(size_t l=0; l<SPECTRUM_SAMPLES; l++){
                            temp[idx1 * SPECTRUM_SAMPLES + l] += Li[l] * (1-remainder1);
                            temp[idx11 * SPECTRUM_SAMPLES + l] += Li[l] * remainder1;
                        }
                    }
                    if(idx2 >= 0){
                        for(size_t l=0; l<SPECTRUM_SAMPLES; l++){
                            temp[(m_M + idx2) * SPECTRUM_SAMPLES + l] += Li[l] * (1-remainder2);
                            temp[(m_M + idx22) * SPECTRUM_SAMPLES + l] += Li[l] * remainder2;
                        }
                    }
                }


                temp[2 * m_M * SPECTRUM_SAMPLES] = rRec.alpha;
                temp[2 * m_M * SPECTRUM_SAMPLES + 1] = 1.0f;
                block->put(samplePos, temp);
                sampler->advance();
            }
        }
    }
    MTS_DECLARE_CLASS()
private:
    Float m_max_distance;
    Float m_min_distance;
};

MTS_IMPLEMENT_CLASS_S(MIFMCWPSDPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIFMCWPSDPathTracer, "MI path tracer");
MTS_NAMESPACE_END
