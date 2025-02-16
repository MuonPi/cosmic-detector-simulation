#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "algebra_types.h"
#include "geometry_types.h"
#include "matrix.h"

constexpr double DEFAULT_EPSILON { 1e-9 };
constexpr double PI { 3.141592653589793238463 };
constexpr double pi() { return PI; }
//template<typename T>
//T const pi = std::acos(-T(1));
constexpr double twopi() { return pi() * 2; }
const double c_sqrt2 { std::sqrt(2.) };

/** @brief calculate euclidean norm (length) of vector
 * @param vec a vector of arbitrary dimension
 * @return the euclidean norm of the vector
 */
auto norm(const Vector& vec) -> double;

/** @brief calculate the cross product of two vectors
 * @param a first vector
 * @param b second vector
 * @return vector containing the cross product of vectors a and b
 * @note the vectors a and b must have dimension=3
 */
auto cross_product(const Vector& a, const Vector& b) -> Vector;

/** @brief check, wether a value is inside an epsilon environment
 * @param value the value to test
 * @param eps optional epsilon value
 * @return true, if the value is within the epsilon range (-eps < value < eps)
 */
bool inEpsilon(double value, double eps = DEFAULT_EPSILON);

/** @brief check, wether two vectors are identical within a epsilon difference
 * @param a the first vector
 * @param b the second vector
 * @param eps optional epsilon value
 * @return true, if the vectors are the same up to a difference of not more than epsilon in its norm
 */
bool isFuzzySame(const std::valarray<double>& a, const std::valarray<double>& b, double eps = DEFAULT_EPSILON);

/** @brief rotate a given point by given angle about given axis
 * @param p the point to be rotated
 * @param rot_axis a vector indicating the rotation axis
 * @param rot_angle the rotation angle in radian
 * @return the point after the rotation operation
 */
Point rotate(const Point& p, const Vector& rot_axis, double rot_angle);

/** @brief get the volume content of a given bounding box
 * @param bb the bounding box given as pair of coordinates (min and max boundaries)
 * @return the volume content
 */
double getBoundingBoxVolume(const std::pair<Point, Point>& bb);

Vector operator*(const matrix2d<double>& lhs, const Vector& rhs);
std::ostream& operator<<(std::ostream& os, const std::valarray<double>& p);
std::ostream& operator<<(std::ostream& os, const matrix2d<double>& m);
bool operator==(const Vector& lhs, const Vector& rhs);
bool operator!=(const Vector& lhs, const Vector& rhs);
