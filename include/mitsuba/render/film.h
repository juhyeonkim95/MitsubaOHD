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

#pragma once
#if !defined(__MITSUBA_RENDER_FILM_H_)
#define __MITSUBA_RENDER_FILM_H_

#include <mitsuba/render/sampler.h>
#include <mitsuba/render/imageblock.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract film base class - used to store samples
 * generated by \ref Integrator implementations.
 *
 * To avoid lock-related bottlenecks when rendering with many cores,
 * rendering threads first store results in an "image block", which
 * is then committed to the film.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Film : public ConfigurableObject {
public:
    /// Ignoring the crop window, return the resolution of the underlying sensor
    inline const Vector2i &getSize() const { return m_size; }

    /// Return the size of the crop window
    inline const Vector2i &getCropSize() const { return m_cropSize; }

    /// Return the offset of the crop window
    inline const Point2i &getCropOffset() const { return m_cropOffset; }

    inline size_t getFrames() const {return m_frames; }

    /// Clear the film
    virtual void clear() = 0;

    /// Merge an image block into the film
    virtual void put(const ImageBlock *block) = 0;

    /// Overwrite the film with the given bitmap and optionally multiply it by a scalar
    virtual void setBitmap(const Bitmap *bitmap, Float multiplier = 1.0f) = 0;

    /// Accumulate a bitmap on top of the radiance values stored in the film
    virtual void addBitmap(const Bitmap *bitmap, Float multiplier = 1.0f) = 0;

    /// Set the target filename (with or without extension)
    virtual void setDestinationFile(const fs::path &filename, uint32_t blockSize) = 0;

    /// Develop the film and write the result to the previously specified filename
    virtual void develop(const Scene *scene, Float renderTime) = 0;

    /**
     * \brief Develop the contents of a subregion of the film and store
     * it inside the given bitmap
     *
     * This may fail when the film does not have an explicit representation
     * of the bitmap in question (e.g. when it is writing to a tiled EXR image)
     *
     * \return \c true upon success
     */
    virtual bool develop(
        const Point2i &offset,
        const Vector2i &size,
        const Point2i &targetOffset,
        Bitmap *target) const = 0;

    /// Does the destination file already exist?
    virtual bool destinationExists(const fs::path &basename) const = 0;

    /**
     * Should regions slightly outside the image plane be sampled to improve
     * the quality of the reconstruction at the edges? This only makes
     * sense when reconstruction filters other than the box filter are used.
     */
    inline bool hasHighQualityEdges() const { return m_highQualityEdges; }

    /// Return whether or not this film records the alpha channel
    virtual bool hasAlpha() const = 0;

    /// Return the image reconstruction filter
    inline ReconstructionFilter *getReconstructionFilter() { return m_filter.get(); }

    /// Return the image reconstruction filter (const version)
    inline const ReconstructionFilter *getReconstructionFilter() const { return m_filter.get(); }

    // =============================================================
    //! @{ \name ConfigurableObject interface
    // =============================================================

    /// Add a child node
    virtual void addChild(const std::string &name, ConfigurableObject *child);

    /// Add an unnamed child
    inline void addChild(ConfigurableObject *child) { addChild("", child); }

    /// Configure the film
    virtual void configure();

    /// Serialize this film to a binary data stream
    virtual void serialize(Stream *stream, InstanceManager *manager) const;

    //! @}
    // =============================================================

    MTS_DECLARE_CLASS()
protected:
    /// Create a film
    Film(const Properties &props);

    /// Unserialize a film
    Film(Stream *stream, InstanceManager *manager);

    /// Virtual destructor
    virtual ~Film();
protected:
    Point2i m_cropOffset;
    Vector2i m_size, m_cropSize;
    bool m_highQualityEdges;
    ref<ReconstructionFilter> m_filter;
    size_t m_frames;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_FILM_H_ */
