#pragma once
#ifndef _FLUIDCANVAS_H_

/**************************************************************************/
#include <cstdint>
#include <opencv2/highgui/highgui.hpp>
#include "Field.h"

/**************************************************************************/
// @name    FluidCanvas
// @brief   Holds every piece of information about the fluid (velocity, density, etc.).
class FluidCanvas {
public:
    // Constructors
    FluidCanvas(
        const std::uint32_t pWidth,
        const std::uint32_t pHeight,
        const double pTimestep,
        const double pDiffusion,
        const double pViscosity);
    ~FluidCanvas();
    
    // Accessors
    std::uint32_t width() const;
    std::uint32_t height() const;
    std::uint32_t total() const;
    
    // Computation
    void step();
    
private:
    // Fluid parameters
    std::uint32_t mWidth;
    std::uint32_t mHeight;
    double mTimestep;
    double mDiffusion;               // Diffusion speed of the fluid
    double mViscosity;               // Fluid viscosity
    
    // Fluid fields
    VectorField mVelocities;         // Holds the velocity at every point of the canvas
    VectorField mBufferField;        // Used to store temporary fields
    
    // Computation functions
    void addDensitySources(VectorField& pVelocities, const VectorField& pSources) const;
};

#endif /* _FLUIDCANVAS_H_ */
