#include "Field.h"

/**************************************************************************/
template <class T>
Field<T>::Field(const std::uint32_t pWidth, const std::uint32_t pHeight, const T& pInitial):
    mWidth{ pWidth },
    mHeight{ pHeight },
    mValues( index(mWidth-1, mHeight-1) + 1, pInitial )
{
}

/**************************************************************************/
template <class T>
Field<T>::Field(const Field<T>& pOther):
    mWidth(pOther.mWidth),
    mHeight(pOther.mHeight),
    mValues(pOther.mValues)
{
}

/**************************************************************************/
template <class T>
Field<T>::Field(Field<T>&& pOther) noexcept:
    mWidth(pOther.mWidth),
    mHeight(pOther.mHeight),
    mValues(std::move(pOther.mValues))
{
    pOther.mWidth = 0;
    pOther.mHeight = 0;
}

/**************************************************************************/
template <class T>
Field<T>::~Field()
{
}

/**************************************************************************/
template <class T>
std::size_t Field<T>::index(const std::uint32_t pX, const std::uint32_t pY) const
{
    return (pY * mWidth + pX);
}

/**************************************************************************/
template <class T>
Field<T>& Field<T>::operator=(const Field<T>& val)
{
    if (this != &val)
    {
        mWidth = val.mWidth;
        mHeight = val.mHeight;
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
        mWidth = val.mWidth;
        mHeight = val.mHeight;
        mValues = std::move(val.mValues);
        
        val.mWidth = 0;
        val.mHeight = 0;
    }
    
    return *this;
}

/**************************************************************************/
template <class T>
Field<T> Field<T>::operator+(const Field<T>& val) const
{
    // TODO
    return *this;
}

/**************************************************************************/
template <class T>
Field<T> Field<T>::operator-(const Field<T>& val) const
{
    // TODO
    return *this;
}

/**************************************************************************/
template <class T>
Field<T> Field<T>::operator*(const double val) const
{
    // TODO
    return *this;
}
