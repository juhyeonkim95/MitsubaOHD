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
#include "../tofpath/path_trace_parts.h"
#include "../tofpath/utils.h"

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
class MIEventExactPathTracer : public MonteCarloIntegrator {
public:
    MIEventExactPathTracer(const Properties &props)
        : MonteCarloIntegrator(props) { 
            m_M = props.getInteger("M", 256);  // number of sampled time
            m_start_time = props.getFloat("start_time", 0.0);
            m_end_time = props.getFloat("end_time", 0.1);
            m_spatial_correlation_mode = props.getString("spatial_correlation_mode", "ray_sampler");
            m_reconnection_threshold = props.getFloat("reconnection_threshold", 0.5);
            m_event_threshold = props.getFloat("event_threshold", 0.1f);
        }

    /// Unserialize from a binary data stream
    MIEventExactPathTracer(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }
    
    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        return Spectrum(0.0);
    }

    std::pair<Spectrum, PathInfo> Li_helper(const RayDifferential &r, RadianceQueryRecord &rRec) const {
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
        PathInfo path_info;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            path_info.add_its(its);
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

            // if (rRec.depth++ >= m_rrDepth) {
            //     /* Russian roulette: try to keep path weights equal to one,
            //        while accounting for the solid angle compression at refractive
            //        index boundaries. Stop with at least some probability to avoid
            //        getting stuck (e.g. due to total internal reflection) */

            //     Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
            //     if (rRec.nextSample1D() >= q)
            //         break;
            //     throughput /= q;
            // }
        }

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;

        return {Li, path_info};
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MIEventExactPathTracer[" << endl
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
        
        sampler->loadSavedState();
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
            auto result = Li_helper(sensorRay, rRec);
            L += result.first * spec;
            sampler->advance();
        }
        L /= N;
        return L;
    }

    std::pair<Float, Float> find_next_event(size_t N, float start_time, 
        float event_threshold, RadianceQueryRecord &rRec, const Sensor *sensor, 
        Sampler *sampler, Point2i& offset) const {
        int count = 0;
        float time = start_time;
        float dt = 0.01;

        Spectrum L0 = single_time(N, start_time, rRec, sensor, sampler, offset);

        float pL = L0[0] * std::exp(event_threshold);
        float a = 1.0;
        float alpha = 1.0;

        // stochastic gradient descent
        while(true){
            float dt = 0.01;

            // newton
            float a_n = a / std::pow(count + 1.0, alpha);// / count;

            // gradient descent
            // float a_n = 0.1 / count;

            //sampler->saveState();
            Spectrum L1 = single_time(4, time, rRec, sensor, sampler, offset);
            //sampler->loadSavedState();
            Spectrum L2 = single_time(4, time + dt, rRec, sensor, sampler, offset);
            
            Float dL = (L2[0] - L1[0]) / dt;
            if(dL == 0.0){
                break;
            }
            
            Float difference = L1[0] - pL;

            // newton
            time = time - 0.5 * difference / dL;
            
            // gradient 
            // time = time - a_n * dL / L0[0] * (difference > 0? 1.0:-1.0);

            if(time < start_time || time > 1.0 || count >= 8){
                break;
            }
            count += 1;
        }
        time = std::max(time, start_time + 1.0f / m_M);
        Spectrum L = single_time(N, time, rRec, sensor, sampler, offset);

        return {time, L[0]};
    }


    std::pair<Float, Float> find_next_event3(size_t N, float start_time, 
        float event_threshold, RadianceQueryRecord &rRec, const Sensor *sensor, 
        Sampler *sampler, Point2i& offset) const {
        int count = 1;
        float time = start_time;

        Spectrum L0 = single_time(N, start_time, rRec, sensor, sampler, offset);
        // Spectrum L1 = single_time(N, start_time, rRec, sensor, sampler, offset);
        
        float pL = L0[0] * std::exp(event_threshold);
        if (L0[0] < 1e-5){
            time = std::max(time, start_time + 1.0f / m_M);
            Spectrum L = single_time(N, time, rRec, sensor, sampler, offset);
            return {time, L[0]};
        }

        while(true){
            float a_n = 1.0 / L0[0] * 1.0 / count;
            Spectrum L = single_time(N, time, rRec, sensor, sampler, offset);

            float new_time = time - a_n * (L[0] - pL);
            new_time = std::max(start_time, std::min(new_time, 1.0f));
            time = new_time;

            count += 1;
            if(count > 1){
                break;
            }
        }
        time = std::max(time, start_time + 1.0f / m_M);
        Spectrum L = single_time(N, time, rRec, sensor, sampler, offset);

        return {time, L[0]};
    }

    std::pair<Float, Float> find_next_event2(size_t N, float start_time, 
        float event_threshold, RadianceQueryRecord &rRec, const Sensor *sensor, 
        Sampler *sampler, Point2i& offset) const {
        float time_epsilon = 0.05;

        Spectrum L0 = single_time(4 * N, start_time, rRec, sensor, sampler, offset);
        Spectrum L1 = single_time(4 * N, start_time + time_epsilon, rRec, sensor, sampler, offset);

        float pL = L0[0] * (1 + event_threshold);

        float t_2 = start_time;
        float t_1 = start_time + time_epsilon;
        float f_2 = L0[0] - pL;
        float f_1 = L1[0] - pL;
        int count = 0;
        while(true){
            if(count >= 5 || f_1 < 1e-10){
                break;
            }
            float t_0 = t_1 - f_1 * (t_1 - t_2) / (f_1 - f_2);
            // if (t_0 < start_time || t_0 > 1.0){
            //     t_0 = 
            // }
            t_0 = std::max(start_time, std::min(t_0, 1.0f));

            float f_0 = single_time(N, t_0, rRec, sensor, sampler, offset)[0] - pL;
            t_2 = t_1;
            t_1 = t_0;
            f_2 = f_1;
            f_1 = f_0;
            count += 1;
        }
        return {t_1, f_1 + pL};
    }

    
    // Overloaded
    void renderBlock(const Scene *scene,
        const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {
        RadianceQueryRecord rRec(scene, sampler);

        block->clear();

        // Float *temp = (Float *) alloca(sizeof(Float) * (m_M * SPECTRUM_SAMPLES + 2));
        size_t N = sampler->getSampleCount();
        Float *temp = (Float *) alloca(sizeof(Float) * (m_M * SPECTRUM_SAMPLES + 2));

        for (size_t i = 0; i<points.size(); ++i) {
            Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
            
            Vector2 jitter = Vector2(0.5);
            Point2 samplePos = Point2(offset) + jitter;
            if (stop)
                break;
            sampler->generate(offset);

            float start_time = 0.0;
            std::vector<float> event_times;
            std::vector<float> event_values;
            sampler->saveState();
            Spectrum L0 = single_time(N, start_time, rRec, sensor, sampler, offset);
            event_times.push_back(0.0);
            event_values.push_back(L0[0]);

            while(true){
                auto event_pos = find_next_event(N, start_time, m_event_threshold, rRec, sensor, sampler, offset);
                auto event_neg = find_next_event(N, start_time, -m_event_threshold, rRec, sensor, sampler, offset);
                
                bool valid_1 = (event_pos.first > start_time && event_pos.first < 1.0);
                bool valid_2 = (event_neg.first > start_time && event_neg.first < 1.0);
                if(!valid_1 && !valid_2){
                    break;
                }
                auto event = event_pos.first < event_neg.first && valid_1 ? event_pos : event_neg;
                start_time = event.first;
                event_times.push_back(event.first);
                event_values.push_back(event.second);
            }

            // if(event_times.size() < 30){
            //     printf("Event Found Size: %d\n", event_times.size());
            //     for(int k=0; k<event_times.size(); k++){
            //         printf("%dth event: %.5f, %.5f\n", k, event_times.at(k), event_values.at(k));
            //     }
            // }
            int j=0;
            int ne = event_values.size();
            for(int m=0; m<m_M; m++){
                float time = m / (float) m_M;
                while(true){
                    if(j+1 >=ne){
                        break;
                    }
                    if (j+1 < ne && time < event_times.at(j+1)){
                        break;
                    }
                    j+=1;
                }
                // while((j <= event_times.size() - 2) && (time > event_values.at(j+1))){
                //     j += 1;
                // }
                Float value;
                if(j == event_times.size() - 1){
                    value = event_values.at(j);
                }
                else{                        
                    Float t1 = event_times.at(j);
                    Float t2 = event_times.at(j + 1);
                    Float f1 = event_values.at(j);
                    Float f2 = event_values.at(j + 1);

                    Float r = (time - t1) / (t2 - t1);
                    value = r * (f2 - f1) + f1;
                }
                temp[m] = value;
            }
            temp[m_M * SPECTRUM_SAMPLES] = rRec.alpha;
            temp[m_M * SPECTRUM_SAMPLES + 1] = 1.0f;
            block->put(samplePos, temp);
        }
    }
    MTS_DECLARE_CLASS()
    
private:
    size_t m_M;
    Float m_start_time;
    Float m_end_time;
    Float m_event_threshold;
    std::string m_spatial_correlation_mode;
    Float m_reconnection_threshold;
};

MTS_IMPLEMENT_CLASS_S(MIEventExactPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MIEventExactPathTracer, "MI path tracer");
MTS_NAMESPACE_END
