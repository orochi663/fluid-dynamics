#pragma once
#ifndef _FIELD_H_
#define _FIELD_H_

/**************************************************************************/
#include <cstdint>
#include <vector>
#include <opencv2/highgui/highgui.hpp>

/**************************************************************************/
// @name    Field
// @brief   Defines a 3D field of values
template <class T>
class Field {
public:
    Field(
        const cv::Vec3i& pSize,
        const T& pInitial);
    Field(const Field<T>& pOther);
    Field(Field<T>&& pOther) noexcept;
    ~Field();
    
    Field<T>& operator=(const Field<T>& val);
    Field<T>& operator=(Field<T>&& val) noexcept;

    Field<T> convolute(const Field<double>& pKernel);

    cv::Vec3i& size() const;
    std::uint32_t total() const;
    
    T& at(const cv::Vec3i& pPosition) const;
    T at(const cv::Vec3i& pPosition);

protected:
    std::size_t index(const cv::Vec3i& pPosition) const;

private:
    cv::Vec3i mSize;
    std::vector<T> mValues;
};

/**************************************************************************/
typedef Field<cv::Vec2d> VectorField;
typedef Field<double> FloatField;

/**************************************************************************/
#include "Field.inc.h"

#endif /* _FIELD_H_ */
