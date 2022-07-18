#pragma once

#include <iostream>
#include <valarray>

#include "matrix.h"

typedef std::valarray<double> Point;
typedef std::valarray<double> Vector;

namespace R3 {
namespace Base {
    const Vector X { 1., 0., 0. };
    const Vector Y { 0., 1., 0. };
    const Vector Z { 0., 0., 1. };
}
const matrix2d<double> Identity { 3,
    { 1., 0., 0.,
        0., 1., 0.,
        0., 0., 1. } };
const Vector NullVec { 0., 0., 0. };
}
