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
        const cv::Vec3u& pSize,
        const T& pInitial);
    Field(const Field<T>& pOther);
    Field(Field<T>&& pOther) noexcept;
    ~Field();
    
    Field<T>& operator=(const Field<T>& val);
    Field<T>& operator=(Field<T>&& val) noexcept;

    Field<T> convolute(const FloatField& pKernel);

    cv::Vec3u& size() const;
    std::uint32_t total() const;
    
    T& at(const cv::Vec3u& pPosition) const inline;
    T at(const cv::Vec3u& pPosition) inline;

protected:
    std::size_t index(const cv::Vec3u& pPosition) const inline;

private:
    std::Vec3u mSize;
    std::vector<T> mValues;
};

/**************************************************************************/
typedef Field<cv::Vec2d> VectorField;
typedef Field<double> FloatField;

/**************************************************************************/
#include "Field.inc.h"

#endif /* _FIELD_H_ */
