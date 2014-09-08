#include "Field.h"
#include <stdexcept>

/**************************************************************************/
template <class T>
Field<T>::Field(const cv::Vec3i& pSize, const T& pInitial):
    mSize{ pSize },
    mValues( total(), pInitial )
{
}

/**************************************************************************/
template <class T>
Field<T>::Field(const Field<T>& pOther):
    mSize(pOther.mSize),
    mValues(pOther.mValues)
{
}

/**************************************************************************/
template <class T>
Field<T>::Field(Field<T>&& pOther) noexcept:
    mSize(pOther.mSize),
    mValues(std::move(pOther.mValues))
{
    pOther.mSize = cv::Vec3i(0, 0, 0);
}

/**************************************************************************/
template <class T>
Field<T>::~Field()
{
}

/**************************************************************************/
template <class T>
std::size_t Field<T>::index(const cv::Vec3i& pPosition) const
{
    if (((pPosition[0] < 0) || (pPosition[0] >= mSize[0])) ||
        ((pPosition[1] < 0) || (pPosition[1] >= mSize[1])) ||
        ((pPosition[2] < 0) || (pPosition[2] >= mSize[2])))
    {
        throw std::out_of_range("Tried to access a cell outside of the field.");
    }

    return (pPosition[2] * mSize[0] * mSize[1]) + 
           (pPosition[1] * mSize[0]) +
            pPosition[0];
}

/**************************************************************************/
template <class T>
Field<T>& Field<T>::operator=(const Field<T>& val)
{
    if (this != &val)
    {
        mSize = val.mSize;
        mValues = val.mValues;
    }
    
    return *this;
}

/**************************************************************************/
template <class T>
Field<T>& Field<T>::operator=(Field<T>&& val) noexcept
{
    if (this != &val)
    {
        mSize = val.mSize;
        mValues = std::move(val.mValues);
        
        val.mSize = cv::Vec3i(0, 0, 0);
    }
    
    return *this;
}

/**************************************************************************/
template <class T>
Field<T> Field<T>::convolute(const Field<double>& pKernel)
{
    Field<T> lResult(mSize, static_cast<T>(0));

    const int tW = mSize[0];
    const int tH = mSize[1];
    const int tD = mSize[2];

    const cv::Vec3i& lKernelSize = pKernel.size();
    
    const uint32_t kW = lKernelSize[0];
    const uint32_t kH = lKernelSize[1];
    const uint32_t kD = lKernelSize[2];
    
    if ((kW%2 == 0) || (kH%2 == 0) || (kD%2 == 0))
        throw std::logic_error("Kernel size should be odd in all directions.");
        
    const int mW = -(kW-1)/2;
    const int mH = -(kH-1)/2;
    const int mD = -(kD-1)/2;
    
    const int MW = mW + kW;
    const int MH = mH + kH;
    const int MD = mD + kD;

    for(int X = 0; X < mSize[0]; ++X)
    {
        for(int Y = 0; Y < mSize[1]; ++Y)
        {
            for(int Z = 0; Z < mSize[2]; ++Z)
            {
                cv::Vec3i lPosition;
                T lValue = static_cast<T>(0);
                
                for(int W = mW; W < MW; ++W)
                {
                    lPosition[0] = X+W;
                    if ((lPosition[0] < 0) || (lPosition[0] >= tW))
                        continue;
                        
                    for(int H = mH; H < MH; ++H)
                    {
                        lPosition[1] = Y+H;
                        if ((lPosition[1] < 0) || (lPosition[1] >= tH))
                            continue;
                        
                        for(int D = mD; D < MD; ++D)
                        {
                            lPosition[2] = Z+D;
                            if ((lPosition[2] < 0) || (lPosition[2] >= tD))
                                continue;
                                
                            // ouf
                            lValue += at(lPosition) * pKernel.at(cv::Vec3i(W, H, D));
                        }
                    }
                }
                
                lPosition = cv::Vec3i(X, Y, Z);
                lResult.at(lPosition) = lValue;
            }
        }
    }
    
    return lResult;
}

/**************************************************************************/
template <class T>
cv::Vec3i& Field<T>::size() const
{
    return mSize;
}

/**************************************************************************/
template <class T>
T& Field<T>::at(const cv::Vec3i& pPosition) const
{
    return mValues.at(index(pPosition));
}

/**************************************************************************/
template <class T>
T Field<T>::at(const cv::Vec3i& pPosition)
{
    return mValues.at(index(pPosition));
}

/**************************************************************************/
template <class T>
std::uint32_t Field<T>::total() const
{
    return mSize[0] * mSize[1] * mSize[2];
}
