#pragma once

#include <cmath>
#include <functional>
#include <valarray>
#include <vector>
#include <stdexcept>

#include "algebra_types.h"

/** @brief Line - struct representing a geometric line in R^n
 * The line is comprised of a reference point p and a direction vector q 
*/
struct Line {
    Point p {};
    Vector q {};
    
    /** @brief evaluate the locus of the line at parameter value t
    * @param t scalar line parameter. t=0 for ref point
    * @return the locus of the line at the selected t-parameter
    */
    auto operator()(double t) const -> Point;

    /** @brief calculate the distance of a given point to the line
    * @param point the point which distance to the line shall be calculated
    * @return norm-distance of the point to the line, i.e. unsigned distance
    */
    auto distance(const Point& point) const -> double;

    /** @brief generate a line object from given ref point and theta, phi angles
    * @param p0 reference point
    * @param theta zenit angle
    * @param phi azimuthal angle
    * @return line object
    */
    static auto generate(Point p0, double theta, double phi) -> Line;
    
    /** @brief rotate the line about given axis by given angle
    * @param rot_axis the axis about which the rotation is performed
    * @param rot_angle the rotation angle (in radians)
    */
    void rotate(const Vector& rot_axis, double rot_angle);
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
    void rotate(const Vector& rot_axis, double rot_angle);
};

class ExtrudedObject {
public:
    ExtrudedObject() = default;
    ExtrudedObject(const std::vector<Point>& vertices, const Point& position, double thickness);
    ExtrudedObject(const Point& position, double radius, double thickness, std::size_t nr_vertices = 32);
    const auto position() const -> Point;
    auto thickness() const -> double;
    auto contains(const Point& point) const -> bool;
    auto intersection(const Line& path) const -> LineSegment;
    auto bounding_box() const -> std::pair<Point, Point>;

private:
    std::vector<Point> m_vertices {};
    Point m_position { 0., 0., 0. };
    double m_thickness { 0. };
    std::vector<Plane> m_planes {};
    matrix2d<double> m_rotation_matrix { R3::Identity };
    auto getPlanes() const -> std::vector<Plane>;
};
