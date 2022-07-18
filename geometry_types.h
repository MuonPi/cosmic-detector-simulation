#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>
#include <valarray>
#include <vector>

#include "algebra_types.h"

/** @brief Line - struct representing a geometric line in R^n
 * The line is comprised of a reference point p and a direction vector q.
 * Any point of the line is represented by the locus p(t) = p + t*q
 * with the line parameter t.
 */
struct Line {
    Point p {}; //<! the reference point of the line
    Vector q {}; //<! the direction vector

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

/** @brief LineSegment - struct representing a geometric line segment in R^n
 * with start and end points indicated by t_start and t_end parameter
 */
struct LineSegment {
    Line line {}; //<! the line defining the segment
    double t_start {}; //<! the line parameter t for the start point of the segment
    double t_end {}; //<! the line parameter t for the end point of the segment
    auto length() const -> double;
};

/** @brief Plane - struct representing a geometric plane in R^n.
 * The plane is comprised of a reference point p and a normal vector.
 */
struct Plane {
    Point p {}; //<! the reference point
    Vector normal {}; //<! the normal vector
    struct no_normal : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    struct no_intersection : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    auto distance(const Point& point) const -> double;
    auto intersection(const Line& line) const -> Point;

    /** @brief Rotate the plane about given axis by given angle.
     * @param rot_axis the axis about which the rotation is performed
     * @param rot_angle the rotation angle (in radians)
     * @note Rotates the normal vector of the plane. The reference point is preserved.
    */
    void rotate(const Vector& rot_axis, double rot_angle);

    /** @brief Rotate the plane through the given rotation matrix.
     * @param rot_matrix the matrix performing the rotation
    */
    void rotate(const matrix2d<double>& rot_matrix);
};

/** @brief ExtrudedObject - class representing a geometric extruded form in R3.
 * The object is defined by a series of 2-dimensional vertex points as outline and a thickness.
 * A global position (the detector's reference point) and a rotation in the local coordinate system around this reference point can be set.
 */
class ExtrudedObject {
public:
    ExtrudedObject() = default;
    ExtrudedObject(const std::vector<Point>& vertices, const Point& position, double thickness);
    ExtrudedObject(const Point& position, double radius, double thickness, std::size_t nr_vertices = 32);
    const auto position() const -> Point;
    void set_position(const Point& new_pos);
    auto thickness() const -> double;
    auto contains(const Point& point) const -> bool;
    auto intersection(const Line& path) const -> LineSegment;
    auto bounding_box() const -> std::pair<Point, Point>;
    void add_rotation(const Vector& rot_axis, double rot_angle);
    auto get_rotation_matrix() -> const matrix2d<double>&;
    void reset_rotation_matrix();

private:
    std::vector<Point> m_vertices {}; //<! the vector of 2d vertices defining the object's outline
    Point m_position { 0., 0., 0. }; //<! the global position and reference point of the object
    double m_thickness { 0. }; //<! the thickness (extrusion) of the object
    std::vector<Plane> m_planes {}; //<! a vector of all surface planes of the object
    matrix2d<double> m_rotation_matrix { R3::Identity }; //<! the 3x3 rotation matrix
    auto getPlanes() const -> std::vector<Plane>;
};
