#include "FluidCanvas.h"

/**************************************************************************/
FluidCanvas::FluidCanvas(const std::uint32_t pWidth, const std::uint32_t pHeight, const double pTimestep, const double pDiffusion, const double pViscosity):
    mWidth(pWidth),
    mHeight(pHeight),
    mTimestep(pTimestep),
    mDiffusion(pDiffusion),
    mViscosity(pViscosity),
    mVelocities(pWidth, pHeight, cv::Vec2d(0.0, 0.0)),
    mBufferField(pWidth, pHeight, cv::Vec2d(0.0, 0.0))
{
}

/**************************************************************************/
FluidCanvas::~FluidCanvas()
{
}

/**************************************************************************/
void FluidCanvas::step()
{

}

/**************************************************************************/
void FluidCanvas::addDensitySources(VectorField& pVelocities, const VectorField& pSources) const
{

}

/**************************************************************************/
std::uint32_t FluidCanvas::width() const
{
    return mWidth;
}

/**************************************************************************/
std::uint32_t FluidCanvas::height() const
{
    return mHeight;
}

/**************************************************************************/
std::uint32_t FluidCanvas::total() const
{
    return mWidth * mHeight;
}
