#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "algebra_types.h"
#include "geometry_types.h"

constexpr double DEFAULT_EPSILON { 1e-9 };
constexpr double pi() { return std::acos(-1); }
constexpr double sqrt2 { std::sqrt(2.) };

std::ostream& operator<<(std::ostream& os, const std::valarray<double>& p);
auto norm(const Vector& vec) -> double;
auto cross_product(const Vector& a, const Vector& b) -> Vector;
bool inEpsilon(double value, double eps = DEFAULT_EPSILON);
bool isFuzzySame(const std::valarray<double>& a, const std::valarray<double>& b);
