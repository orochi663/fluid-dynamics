#include "FluidCanvas.h"

/**************************************************************************/
FluidCanvas::FluidCanvas(const cv::Vec3i& pSize, const double pTimestep, const double pDiffusion, const double pViscosity):
    mSize(pSize),
    mTimestep(pTimestep),
    mDiffusion(pDiffusion),
    mViscosity(pViscosity),
    mVelocities(pSize, cv::Vec2d(0.0, 0.0)),
    mBufferField(pSize, cv::Vec2d(0.0, 0.0))
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
const cv::Vec3i& FluidCanvas::size() const
{
    return mSize;
}

/**************************************************************************/
std::uint32_t FluidCanvas::total() const
{
    return mSize[0] * mSize[1] * mSize[2];
}
