#pragma once
#ifndef _FLUIDCANVAS_H_

/**************************************************************************/
#include <cstdint>
#include <opencv2/highgui/highgui.hpp>
#include "VectorField.h"

/**************************************************************************/
// @name    FluidCanvas
// @brief   Holds every piece of information about the fluid (velocity, density, etc.).
class FluidCanvas {
public:
    typedef std::uint32_t Length;

    FluidCanvas(
        const Length pWidth,
        const Length pHeight,
        const double pTimestep,
        const double pDiffusion,
        const double pViscosity);
    ~FluidCanvas();
    
private:
    // Fluid parameters
    Length mWidth;
    Length mHeight;
    double mTimestep;
    double mDiffusion;          // Diffusion speed of the fluid
    double mViscosity;          // Fluid viscosity
    
    // Fluid fields
    VectorField mVelocities;    // Holds the velocity at every point of the canvas
    VectorField mBufferField;   // Used to store temporary fields
};

#endif /* _FLUIDCANVAS_H_ */
