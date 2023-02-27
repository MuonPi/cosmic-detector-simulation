#pragma once

#include "algebra_utils.h"
#include <cmath>
#include <vector>

template <typename T>
struct DataItem {
    T value {};
    T error {};
};

template <typename T, typename U>
using MeasurementVector = std::vector<std::pair<DataItem<T>, DataItem<U>>>;

constexpr double toDeg(double x) { return x * 180 / pi(); }
constexpr double toRad(double x) { return x * pi() / 180; }

void export_file(const MeasurementVector<double, double>& data, const std::string& filename);
