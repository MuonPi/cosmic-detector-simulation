#pragma once

#include <cmath>
#include <functional>
#include <valarray>
#include <vector>

#include "algebra_types.h"

struct Line {
    Point p{};
    Vector q{};
    auto operator()(double t) const -> Point;
    auto distance(const Point &point) const -> double;
    static auto generate(Point p0, double theta, double phi) -> Line;
};

struct LineSegment {
    Line line {};
    double t_start {};
    double t_end {};
    auto length() const -> double;
};

struct Plane {
    Point p{};
    Vector q{};
    Vector r{};
    auto operator()(double t, double s) const -> Point;
    auto normal() const -> Vector;
    auto distance(const Point &point) const -> double;
    auto intersection(const Line &line) const -> Point;
};

class ExtrudedObject {
public:
    ExtrudedObject() = default;
    ExtrudedObject(const std::vector<Point> &vertices, const Point& position, double thickness);
    auto contains(const Point &point) const -> bool;
    auto intersection(const Line &path) const -> LineSegment;
    auto bounding_box() const -> std::pair<Point,Point>;
private:
    std::vector<Point> m_vertices{};
    Point m_position {0., 0., 0.};
    double m_thickness{0.};
    std::vector<Plane> m_planes {};
  
    auto getPlanes() const -> std::vector<Plane>;
};

