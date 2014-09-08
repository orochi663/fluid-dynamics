#pragma once
#ifndef _FIELD_H_
#define _FIELD_H_

/**************************************************************************/
#include <cstdint>
#include <vector>
#include <opencv2/highgui/highgui.hpp>

/**************************************************************************/
// @name    Field
// @brief   Defines a 2D field of values
template <class T>
class Field {
public:
    Field(
        const std::uint32_t pWidth,
        const std::uint32_t pHeight,
        const T& pInitial);
    Field(const Field<T>& pOther);
    Field(Field<T>&& pOther) noexcept;
    ~Field();
    
    Field<T>& operator=(const Field<T>& val);
    Field<T>& operator=(Field<T>&& val) noexcept;
    
    Field<T> operator+(const Field<T>& val) const;
    Field<T> operator-(const Field<T>& val) const;
    Field<T> operator*(const double val) const;

protected:
    std::size_t index(
        const std::uint32_t pX,
        const std::uint32_t pY) const;

private:
    std::uint32_t mWidth;
    std::uint32_t mHeight;
    std::vector<T> mValues;
};

/**************************************************************************/
typedef Field<cv::Vec2d> VectorField;
typedef Field<double> FloatField;

/**************************************************************************/
#include "Field.inc.h"

#endif /* _FIELD_H_ */
