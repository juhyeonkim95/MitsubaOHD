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
#include "event_utils.h"
#include "../tofpath/path_trace_parts.h"

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Event Path tracer", "Average path length", EAverage);

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
class MIEventCorrelatedRoughPathTracer : public MonteCarloIntegrator {
public:
    MIEventCorrelatedRoughPathTracer(const Properties &props)
        : MonteCarloIntegrator(props) { 
            m_M = props.getInteger("M", 256);  // number of sampled time
            m_start_time = props.getFloat("start_time", 0.0);
            m_end_time = props.getFloat("end_time", 0.1);
            m_event_threshold = props.getFloat("event_threshold", 0.1f);
        }

    /// Unserialize from a binary data stream
    MIEventCorrelatedRoughPathTracer(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);

            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

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

            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                (bsdf->getType() & BSDF::ESmooth)) {
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
                        Li += throughput * value * bsdfVal * weight;
                    }
                }
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
            eta *= bRec.eta;

            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                    scene->pdfEmitterDirect(dRec) : 0;
                Li += throughput * value * miWeight(bsdfPdf, lumPdf);
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

        return Li;
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
        oss << "MIEventCorrelatedRoughPathTracer[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    Spectrum single_time(size_t N, float time_ratio, RadianceQueryRecord &rRec, const Sensor *sensor, Sampler *sampler, Point2i& offset) const{
        Spectrum L(0.0f);

        Float diffScaleFactor = 1.0f /
            std::sqrt((Float) sampler->getSampleCount());
        
        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;
        uint32_t queryType = RadianceQueryRecord::ESensorRay;
        if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
            queryType &= ~RadianceQueryRecord::EOpacity;

        bool needsApertureSample = sensor->needsApertureSample();
        Float time = time_ratio * (m_end_time - m_start_time) + m_start_time;
        
        for (size_t j = 0; j<N; j++) {
            Vector2 jitter = Vector2(rRec.nextSample2D()) * 2.0 - Vector2(1.0);
            Point2 samplePos = Point2(offset) + jitter;
            if (needsApertureSample)
                apertureSample = rRec.nextSample2D();
            rRec.newQuery(queryType, sensor->getMedium());
            Spectrum spec = sensor->sampleRayDifferential(
                sensorRay, samplePos, apertureSample, time);

            sensorRay.scaleDifferential(diffScaleFactor);
            sensorRay.setTime(time);
            Spectrum result = Li(sensorRay, rRec);
            L += result * spec;
            sampler->advance();
        }
        L /= N;
        return L;
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

        ref<Sampler> sampler2 = sampler->clone();
        RadianceQueryRecord rRec2(scene, sampler2);

        Point2 apertureSample(0.5f);
        Float timeSample = 0.5f;
        RayDifferential sensorRay;

        block->clear();

        uint32_t queryType = RadianceQueryRecord::ESensorRay;

        if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
            queryType &= ~RadianceQueryRecord::EOpacity;

        Float *temp = (Float *) alloca(sizeof(Float) * (m_M * SPECTRUM_SAMPLES + 2));

        
        for (size_t i = 0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            if (stop)
                break;

            sampler->generate(offset);
            // clear buffer
            for(size_t j=0; j<(m_M * SPECTRUM_SAMPLES + 2); j++){
                temp[j] = 0.0;
            }
            sampler->saveState();
            sampler2->saveState();
            int N = sampler->getSampleCount();
            float dt = 1.0 / m_M;

            // (1) Ours
            // (a) build rough distribution
            int M_small = m_M / 8;
            int N_small = N;
            PiecewiseLinearSamples samples;
            for(int i=0; i < M_small; i++){
                sampler->loadSavedState();
                Float time = (i + 0.5) / M_small;
                Spectrum L = single_time(N_small, time, rRec, sensor, sampler, offset);
                samples.add_sample(time, L[0]);
            }

            // (b) search discontinuity
            Intersection its_prev;
            std::vector<float> discont_times;
            for(int i=0; i<m_M; i++){
                Float time_ratio = (i + 0.5) / m_M;

                Float diffScaleFactor = 1.0f /
                    std::sqrt((Float) sampler->getSampleCount());
                
                Point2 apertureSample(0.5f);
                Float timeSample = 0.5f;
                RayDifferential sensorRay;
                uint32_t queryType = RadianceQueryRecord::ESensorRay;
                if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
                    queryType &= ~RadianceQueryRecord::EOpacity;

                bool needsApertureSample = sensor->needsApertureSample();
                Float time = time_ratio * (m_end_time - m_start_time) + m_start_time;
                
                Vector2 jitter = Vector2(0.5f);
                Point2 samplePos = Point2(offset) + jitter;
                if (needsApertureSample)
                    apertureSample = rRec.nextSample2D();
                rRec.newQuery(queryType, sensor->getMedium());
                Spectrum spec = sensor->sampleRayDifferential(
                    sensorRay, samplePos, apertureSample, time);

                sensorRay.scaleDifferential(diffScaleFactor);
                sensorRay.setTime(time);
                rRec.rayIntersect(sensorRay);

                if(i>0)
                {
                    bool consistency = check_its_consistency(its_prev, rRec.its);
                    if(!consistency){
                        discont_times.push_back(time_ratio);
                    }
                }
                its_prev = rRec.its;
            }
            
            PiecewiseLinearSamples new_samples;
            for(int i=0; i<discont_times.size(); i++){
                sampler->loadSavedState();
                Float time = discont_times.at(i);
                Spectrum L = single_time(N_small, time-dt, rRec, sensor, sampler, offset);
                new_samples.add_sample(time-dt, L[0]);
                sampler->loadSavedState();
                L = single_time(N_small, time, rRec, sensor, sampler, offset);
                new_samples.add_sample(time, L[0]);
            }
            samples.merge(new_samples);


            float time = 0.0;
            Spectrum L0 = single_time(N * 16, time, rRec, sensor, sampler, offset);

            // (b) find events
            for(int i=0; i<4; i++){
                PiecewiseLinearSamples new_samples;
                float eventtime = 0.0;
                Float target_L = L0[0];
                while(true){
                    Float eventtime1 = samples.get_time(target_L * std::exp(+m_event_threshold), eventtime);
                    Float eventtime2 = samples.get_time(target_L * std::exp(-m_event_threshold), eventtime);
                    eventtime = std::min(eventtime1, eventtime2);
                    if(eventtime > 1.0){
                        break;
                    }
                    sampler->loadSavedState();
                    Spectrum newL = single_time(N, eventtime, rRec, sensor, sampler, offset);
                    new_samples.add_sample(eventtime, newL[0]);
                    target_L = eventtime1 < eventtime2 ? target_L * std::exp(+m_event_threshold) : target_L * std::exp(-m_event_threshold);
                }
                samples.merge(new_samples);
            }
            samples.put_to(temp, m_M);

            // (1-1) ESIM Variant
            // float time = 0.0;
            // Spectrum L0 = single_time(N * 16, time, rRec, sensor, sampler, offset);
            // PiecewiseLinearSamples samples;
            // float a_n = 0.5;

            // while(true){
            //     Float pL = L0[0] * std::exp(+m_event_threshold);
            //     Float nL = L0[0] * std::exp(-m_event_threshold);

            //     sampler->loadSavedState();
            //     Spectrum L1 = single_time(N, time, rRec, sensor, sampler, offset);
            //     sampler->loadSavedState();
            //     Spectrum L2 = single_time(N, time + dt, rRec, sensor, sampler, offset);
                
            //     Float L = (safe_log(L2[0]) + safe_log(L1[0])) / 2;
            //     Float dL = (safe_log(L2[0]) - safe_log(L1[0])) / dt;
                
            //     Float difference = m_event_threshold;
            //     samples.add_sample(time, (L1[0] + L2[0]) / 2);

            //     if(dL != 0.0){
            //         float new_time = time + a_n * std::abs(difference / dL);
            //         time = std::max(new_time, time + dt);
            //     } else {
            //         time = time + dt;
            //     }
            //     if(time > 1.0){
            //         break;
            //     }
            // }
            // samples.put_to(temp, m_M);

            // (1) ESIM
            // float time = 0.0;
            // PiecewiseLinearSamples samples;

            // while(true){
            //     float dt = 0.01;
            //     sampler->loadSavedState();
            //     Spectrum L1 = single_time(N / 2, time, rRec, sensor, sampler, offset);
            //     sampler->loadSavedState();
            //     Spectrum L2 = single_time(N / 2, time + dt, rRec, sensor, sampler, offset);
                
            //     Float L = (safe_log(L1[0]) + safe_log(L2[0])) / 2;
            //     Float dL = (safe_log(L2[0]) - safe_log(L1[0])) / dt;
                
            //     Float difference = m_event_threshold;
            //     samples.add_sample(time, (L1[0] + L2[0]) / 2);

            //     if(dL != 0.0){
            //         float new_time = time + 0.5 * std::abs(difference / dL);
            //         time = std::max(new_time, time + 1.0f / m_M);
            //     } else {
            //         time = time + 1.0f / m_M;
            //     }
                
            //     if(time > 1.0){
            //         break;
            //     }
            // }
            // samples.put_to(temp, m_M);
            

            // (2) uniform
            // for(size_t m=0; m < m_M; m++){
            //     sampler->loadSavedState();
            //     Spectrum L = single_time(N, (m + 0.5f) / m_M, rRec, sensor, sampler, offset);

            //     for(size_t l=0; l<SPECTRUM_SAMPLES; l++){
            //         temp[m * SPECTRUM_SAMPLES + l] = L[l];
            //     }
            // }

            temp[m_M * SPECTRUM_SAMPLES] = rRec.alpha;
            temp[m_M * SPECTRUM_SAMPLES + 1] = 1.0f;

            Point2 samplePos = Point2(offset);

            block->put(Point2(offset), temp);
        }
    }
    MTS_DECLARE_CLASS()
private:
    size_t m_M;
    Float m_start_time;
    Float m_end_time;
    Float m_event_threshold;
};

MTS_IMPLEMENT_CLASS_S(MIEventCorrelatedRoughPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIEventCorrelatedRoughPathTracer, "MI path tracer");
MTS_NAMESPACE_END
