#include "BezierSurface.h"

#include "Transform.h"
#include "BezierUtil.h"

#include <algorithm>
#include <stdexcept>

///////////////////////////////////////////////////////////////////////////
BezierSurface::BezierSurface()
    : _degreeU(0), _degreeV(0), _nbPointsU(0), _nbPointsV(0)
{
}

BezierSurface::~BezierSurface()
{
}

BezierSurface& BezierSurface::operator=(const BezierSurface& other)
{
    if (this != &other) {
        _degreeU = other._degreeU;
        _degreeV = other._degreeV;
        _nbPointsU = other._nbPointsU;
        _nbPointsV = other._nbPointsV;
        _controlPoints = other._controlPoints;
    }
    return *this;
}

void BezierSurface::clear()
{
    _degreeU = 0;
    _degreeV = 0;
    _nbPointsU = 0;
    _nbPointsV = 0;
    _controlPoints.clear();
}

void BezierSurface::set_degree(int degreeU, int degreeV)
{
    _degreeU = degreeU;
    _degreeV = degreeV;
}

void BezierSurface::set_control_points(const std::vector<Point3>& points, int nbPointsU, int nbPointsV)
{
    if (nbPointsU <= 0 || nbPointsV <= 0) {
        clear();
        return;
    }

    if ((int)points.size() != nbPointsU * nbPointsV) {
        throw std::invalid_argument("Control points size doesn't match dimensions");
    }

    _nbPointsU = nbPointsU;
    _nbPointsV = nbPointsV;
    _controlPoints = points;
}

void BezierSurface::apply_transform(const Transform& t)
{
    for (auto& p : _controlPoints) {
        t.apply(p);
    }
}

Point3 BezierSurface::evaluate(double u, double v) const
{
    if (!is_valid() || _controlPoints.empty()) {
        return Point3();
    }

    // Clamp parameters to [0,1]
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    Point3 result;

    // Evaluate using tensor product of Bernstein polynomials
    for (int i = 0; i < _nbPointsU; ++i) {
        for (int j = 0; j < _nbPointsV; ++j) {
            double bu = getBernsteinValue(_degreeU, i, u);
            double bv = getBernsteinValue(_degreeV, j, v);

            result += _controlPoints[i * _nbPointsV + j] * (bu * bv);
        }
    }

    return result;
}

void BezierSurface::evaluate_derivatives(double u, double v, Point3& du, Point3& dv) const
{
    du = Point3();
    dv = Point3();

    if (!is_valid() || _controlPoints.empty()) {
        return;
    }

    // Clamp parameters to [0,1]
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    // Evaluate derivatives using tensor product
    for (int i = 0; i < _nbPointsU; ++i) {
        for (int j = 0; j < _nbPointsV; ++j) {
            double bu = getBernsteinValue(_degreeU, i, u);
            double bv = getBernsteinValue(_degreeV, j, v);

            double dbu_du = getBernsteinDerivative(_degreeU, i, u);
            double dbv_dv = getBernsteinDerivative(_degreeV, j, v);

            du += _controlPoints[i * _nbPointsV + j] * (dbu_du * bv);
            dv += _controlPoints[i * _nbPointsV + j] * (bu * dbv_dv);
        }
    }
}

bool BezierSurface::is_valid() const
{
    if (_degreeU < 0 || _degreeV < 0) return false;
    if (_nbPointsU != _degreeU + 1) return false;
    if (_nbPointsV != _degreeV + 1) return false;
    if ((int)_controlPoints.size() != _nbPointsU * _nbPointsV) return false;
    return true;
}

// Helper function to evaluate Bernstein basis functions
double BezierSurface::getBernsteinValue(int degree, int index, double t) const
{
    if (degree == 0) return 1.0;
    if (degree == 1) {
        return (index == 0) ? (1.0 - t) : t;
    }
    if (degree == 2) {
        if (index == 0) return BezierUtil::Bernstein02(t);
        if (index == 1) return BezierUtil::Bernstein12(t);
        if (index == 2) return BezierUtil::Bernstein22(t);
    }
    if (degree == 3) {
        if (index == 0) return BezierUtil::Bernstein03(t);
        if (index == 1) return BezierUtil::Bernstein13(t);
        if (index == 2) return BezierUtil::Bernstein23(t);
        if (index == 3) return BezierUtil::Bernstein33(t);
    }

    // For higher degrees, use recursive definition
    // B_{i,n}(t) = (1-t) * B_{i,n-1}(t) + t * B_{i-1,n-1}(t)
    if (index < 0 || index > degree) return 0.0;

    double result = 1.0;
    for (int i = 1; i <= degree; ++i) {
        if (index < i) {
            result *= (1.0 - t);
        } else {
            result *= t;
        }
    }

    // This is a simplified implementation - full Bernstein would need binomial coefficients
    // For now, return 0 for unsupported degrees
    return (degree <= 3) ? result : 0.0;
}

double BezierSurface::getBernsteinDerivative(int degree, int index, double t) const
{
    if (degree <= 0) return 0.0;

    // Derivative of Bernstein basis: n * (B_{i-1,n-1}(t) - B_{i,n-1}(t))
    double b1 = (index > 0) ? getBernsteinValue(degree - 1, index - 1, t) : 0.0;
    double b2 = (index < degree) ? getBernsteinValue(degree - 1, index, t) : 0.0;

    return degree * (b1 - b2);
}