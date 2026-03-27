#ifndef BezierSurface_
#define BezierSurface_

#include <vector>
#include "Geometry.h"

class Transform;

///////////////////////////////////////////////////////////////////////////
class BezierSurface
{
public:
    BezierSurface();
    virtual ~BezierSurface();
    BezierSurface& operator=(const BezierSurface& other);

    void clear();

    void set_degree(int degreeU, int degreeV);
    int degree_u() const { return _degreeU; }
    int degree_v() const { return _degreeV; }

    void set_control_points(const std::vector<Point3>& points, int nbPointsU, int nbPointsV);
    const std::vector<Point3>& control_points() const { return _controlPoints; }
    std::vector<Point3>& control_points() { return _controlPoints; }

    int nb_control_points_u() const { return _nbPointsU; }
    int nb_control_points_v() const { return _nbPointsV; }
    int nb_control_points() const { return _controlPoints.size(); }

    void apply_transform(const Transform& t);

    // Evaluate point on surface
    Point3 evaluate(double u, double v) const;

    // Evaluate derivatives
    void evaluate_derivatives(double u, double v, Point3& du, Point3& dv) const;

    // Check if valid (degrees match control point counts)
    bool is_valid() const;

private:
    double getBernsteinValue(int degree, int index, double t) const;
    double getBernsteinDerivative(int degree, int index, double t) const;

    int _degreeU, _degreeV;
    int _nbPointsU, _nbPointsV;
    std::vector<Point3> _controlPoints;
};
///////////////////////////////////////////////////////////////////////////

#endif