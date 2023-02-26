#pragma once

#include <vector>
#include <cmath>

constexpr double toDeg(double x) { return x * 180 / pi(); }
constexpr double toRad(double x) { return x * pi() / 180; }

/** @brief The cumulative distribution function (CDF) for the cos^2(x) PDF distribution
 * This CDF is used for the calculation of the Probability Density Function (PDF)
 * of the cos^2(x) distribution for generating random values following the angular distribution
 * of the muon tracks
 */
auto cos2cdf = [](double x) {
    //return cdf to following pdf: cos^2(x)
    return (2 / pi()) * (x / 2. + sin(2. * x) / 4.) + 0.5; //from Wolfram Alpha
};

template <typename T>
void Append(std::vector<T>& a, const std::vector<T>& b)
{
    a.reserve(a.size() + b.size());
    a.insert(a.end(), b.begin(), b.end());
}
