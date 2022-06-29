#include <vector>
#include <cmath>

#include "algebra_types.h"
#include "geometry_types.h"
#include "algebra_utils.h"

std::ostream &operator<<(std::ostream &os, const std::valarray<double> &p) {
  os << "(";
  for (double element : p) {
    os << element << " ";
  }
  os << ")";
  return os;
}

auto norm(const Vector &vec) -> double {
  constexpr double exponent{2};
  return std::sqrt(std::pow(vec, exponent).sum());
}

auto cross_product(const Vector &a, const Vector &b) -> Vector {
  //    static_assert(a.size() == 3);
  //    static_assert(b.size() == 3);
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

bool inEpsilon(double value, double eps) {
    if ( std::fabs(value) > eps ) return false;
    return true;
}

bool isFuzzySame(const std::valarray<double>& a, const std::valarray<double>& b) {
    return inEpsilon(norm(a-b));
}
