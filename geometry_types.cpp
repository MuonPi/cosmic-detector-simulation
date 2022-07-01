#include "geometry_types.h"
#include "algebra_utils.h"

#include <iostream>
#include <limits>

auto Line::operator()(double t) const -> Point
{
    return p + t * q;
}

auto Line::distance(const Point& point) const -> double
{
    return norm((point - p) - q * q * (point - p));
}

auto LineSegment::length() const -> double
{
    return norm(line(t_start) - line(t_end));
}

auto Plane::operator()(double t, double s) const -> Point
{
    return p + t * q + s * r;
}

auto Plane::normal() const -> Vector
{
    Vector normalvec { cross_product(q, r) };
    double veclength { norm(normalvec) };
    if (inEpsilon(veclength)) {
        throw no_normal("plane has degenerate normal");
        return { };
        std::cerr << "Error in Plane::normal(): normal vector with length=0!\n";
        return { };
    }
    return normalvec / veclength;
}

auto Plane::distance(const Point& point) const -> double
{
    return ((point - p) * normal()).sum();
}

auto Plane::intersection(const Line& line) const -> Point
{
    if (isFuzzySame(p, line.p)) {
        // reference points are identical
        return p;
    }
    double v1 { ((p - line.p) * normal()).sum() };
    double v2 { (line.q * normal()).sum() };

    if (inEpsilon(v1)) {
        //std::cerr << "Error in Plane::intersection(const Line&): line is contained entirely in plane v1="<<v1<<" v2="<<v2<<" this->p="<<p<<" line.p="<<line.p<<" this->normal()="<<normal()<<" (p-line.p)="<<p-line.p<<" (p-line.p)*normal="<<(p-line.p)*normal()<<"\n";
        return p;
    }

    if (inEpsilon(v2)) {
        //std::cerr << "Error in Plane::intersection(const Line&): line is not intersecting the plane\n";
        throw no_intersection("Plane::intersection(const Line&): line is not intersecting the plane");
        return { };
    }
    double t { v1 / v2 };
    return line(t);
}

ExtrudedObject::ExtrudedObject(const std::vector<Point>& vertices, const Point& position, double thickness)
    : m_vertices(vertices)
    , m_position(position)
    , m_thickness(thickness)
{
    m_planes = getPlanes();
}

ExtrudedObject::ExtrudedObject(const Point& position, double radius, double thickness, std::size_t nr_vertices)
    : m_position(position)
    , m_thickness(thickness)
{
    for ( std::size_t n { 0 }; n < nr_vertices; ++n ) {
        const double angle { twopi() * n / nr_vertices };
        Point p {
            m_position[0] + radius * std::cos(angle),
            m_position[1] + radius * std::sin(angle)
        };
        m_vertices.emplace_back(std::move(p));
    }
    m_planes = getPlanes();
}


auto ExtrudedObject::contains(const Point& point) const -> bool
{
    for (const auto& plane : m_planes) {
        double dist { plane.distance(point) };
        if (dist > DEFAULT_EPSILON)
            return false;
    }
    return true;
}

auto ExtrudedObject::intersection(const Line& path) const -> LineSegment
{
    std::vector<Point> hitpoints {};
    for (const auto& plane : m_planes) {
        Point hitpoint {};
        try {
            hitpoint = { plane.intersection(path) };
        } catch (std::exception& e) {
            continue;
        }
        if (hitpoint[2] < m_position[2] || hitpoint[2] > m_position[2] + m_thickness)
            continue;
        hitpoints.push_back(std::move(hitpoint));
    }
    if (hitpoints.empty()) {
        //std::cerr<<"Error in ExtrudedObject::intersection(const Line&): no intersections with volume found!\n";
        return LineSegment {};
    }
    auto it = hitpoints.begin();
    while (it != hitpoints.end()) {
        if (contains(*it)) {
            //std::cout<<"found intersection point "<<*it<<"\n";
            ++it;
        } else
            it = hitpoints.erase(it);
    }
    if (hitpoints.size() > 1) {
        it = hitpoints.begin();
        while (it != hitpoints.end()) {
            auto it2 = std::next(it);
            while (it2 != hitpoints.end()) {
                if (isFuzzySame(*it, *it2)) {
                    it2 = hitpoints.erase(it2);
                } else
                    ++it2;
            }
            //std::cout<<"found intersection point "<<*it<<"\n";
            ++it;
        }
    }
    if (hitpoints.size() == 2) {
        return LineSegment {
            { hitpoints[0], hitpoints[1] - hitpoints[0] },
            0., 1.
        };
    } else if (hitpoints.size() > 2) {
        //        std::cerr<<"ExtrudedObject::intersection(const Line&): strange nr of intersection points: " << hitpoints.size() << "\n";
    }
    return LineSegment {};
}

auto ExtrudedObject::getPlanes() const -> std::vector<Plane>
{
    std::vector<Plane> planes {};
    if (m_vertices.size() < 3) {
        throw std::runtime_error("Error in ExtrudedObject::getPlanes(): insufficient number of vertices ("+ std::to_string(m_vertices.size())+")!");
        std::cerr << "Error in ExtrudedObject::getPlanes(): insufficient number of vertices (" << m_vertices.size() << ")!\n";
        return planes;
    }
    for (auto vertex { m_vertices.begin() };
         vertex != std::prev(m_vertices.end());
         ++vertex) {
        Point p0 { Point { (*vertex)[0], (*vertex)[1], 0. } + m_position };
        Point p1 { Point { (*std::next(vertex))[0], (*std::next(vertex))[1], 0. } + m_position };
        Point p2 { Point { p0[0], p0[1], m_thickness } + m_position };
        Plane plane { p0, p1 - p0, p2 - p0 };
        planes.push_back(std::move(plane));
    }
    Point p0 { Point { (*std::prev(m_vertices.end()))[0], (*std::prev(m_vertices.end()))[1], 0. } + m_position };
    Point p1 { Point { (*m_vertices.begin())[0], (*m_vertices.begin())[1], 0. } + m_position };
    Point p2 { Point { p0[0], p0[1], m_thickness } + m_position };
    Plane plane { p0, p1 - p0, p2 - p0 };
    planes.push_back(std::move(plane));
    planes.push_back({ m_position, { -1., 0., 0. }, { 0., 1., 0. } });
    planes.push_back({ m_position + Vector { 0., 0., m_thickness }, { 1., 0., 0. }, { 0., 1., 0. } });
    return planes;
}

auto ExtrudedObject::bounding_box() const -> std::pair<Point, Point>
{
    Point min_coordinates {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        0.
    };
    Point max_coordinates {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        0.
    };
    for (const auto& vertex : m_vertices) {
        if (vertex[0] < min_coordinates[0]) {
            min_coordinates[0] = vertex[0];
        }
        if (vertex[1] < min_coordinates[1]) {
            min_coordinates[1] = vertex[1];
        }
        if (vertex[0] > max_coordinates[0]) {
            max_coordinates[0] = vertex[0];
        }
        if (vertex[1] > max_coordinates[1]) {
            max_coordinates[1] = vertex[1];
        }
    }
    min_coordinates += m_position;
    max_coordinates += m_position;
    max_coordinates += { 0., 0., m_thickness };
    return std::make_pair<Point, Point>(std::move(min_coordinates), std::move(max_coordinates));
}

auto Line::generate(Point p0, double theta, double phi) -> Line
{
    return {
        p0,
        { std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta) }
    };
}
