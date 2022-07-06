#pragma once

#include <cmath>
#include <functional>
#include <valarray>
#include <vector>
#include <stdexcept>

#include "algebra_types.h"

struct Line {
    Point p {};
    Vector q {};
    auto operator()(double t) const -> Point;
    auto distance(const Point& point) const -> double;
    static auto generate(Point p0, double theta, double phi) -> Line;
};

struct LineSegment {
    Line line {};
    double t_start {};
    double t_end {};
    auto length() const -> double;
};

struct Plane {
    Point p {};
    Vector q {};
    Vector r {};
    struct no_normal : std::runtime_error { using std::runtime_error::runtime_error; };
    struct no_intersection : std::runtime_error { using std::runtime_error::runtime_error; };
    auto operator()(double t, double s) const -> Point;
    auto normal() const -> Vector;
    static auto fromNormalVector(const Point& ref_point, const Vector& vec) -> Plane;
    auto distance(const Point& point) const -> double;
    auto intersection(const Line& line) const -> Point;
};

class ExtrudedObject {
public:
    ExtrudedObject() = default;
    ExtrudedObject(const std::vector<Point>& vertices, const Point& position, double thickness);
    ExtrudedObject(const Point& position, double radius, double thickness, std::size_t nr_vertices = 32);
    auto contains(const Point& point) const -> bool;
    auto intersection(const Line& path) const -> LineSegment;
    auto bounding_box() const -> std::pair<Point, Point>;

private:
    std::vector<Point> m_vertices {};
    Point m_position { 0., 0., 0. };
    double m_thickness { 0. };
    std::vector<Plane> m_planes {};

    auto getPlanes() const -> std::vector<Plane>;
};
