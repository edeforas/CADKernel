#include "NurbsCurveJoin.h"
#include <cmath>

static bool build_linear_connector(const Point3& p0, const Point3& p3, NurbsCurve& out)
{
	out.clear();
	std::vector<Point3> pts = { p0, p3 };
	out.set_degree(1);
	out.set_points(pts);
	out.set_equals_weights();
	out.set_knots({ 0., 0., 1., 1. });
	return true;
}

static bool build_cubic_transition(const Point3& p0, const Point3& p3, const Point3& t0, const Point3& t1, double h, NurbsCurve& out)
{
	Point3 q1 = p0 + t0 * h;
	Point3 q2 = p3 - t1 * h;
	std::vector<Point3> pts = { p0, q1, q2, p3 };
	out.set_degree(3);
	out.set_points(pts);
	out.set_equals_weights();
	out.set_knots({ 0., 0., 0., 0., 1., 1., 1., 1. });
	return true;
}

static bool build_quartic_fillet(const Point3& p0, const Point3& p3, const Point3& t0, const Point3& t1, double d, double h, NurbsCurve& out)
{
	Point3 q1 = p0 + t0 * (h * 0.5);
	Point3 q3 = p3 - t1 * (h * 0.5);
	Point3 q2 = (p0 + p3) * 0.5 + (t0 - t1) * (d * 0.05);
	std::vector<Point3> pts = { p0, q1, q2, q3, p3 };
	out.set_degree(4);
	out.set_points(pts);
	out.set_equals_weights();
	out.set_knots({ 0.,0.,0.,0.,0.,1.,1.,1.,1.,1. });
	return true;
}

bool NurbsCurveJoin::create_connector(const NurbsCurve& c1, const NurbsCurve& c2, Continuity continuity, NurbsCurve& out)
{
	if (c1.nb_points() < 2 || c2.nb_points() < 2)
		return false;

	Point3 p0 = c1.points().back();
	Point3 p3 = c2.points().front();
	out.clear();

	if (continuity == G0)
		return build_linear_connector(p0, p3, out);

	Point3 t0, t1;
	if (!c1.tangent(1.0, t0) || !c2.tangent(0.0, t1))
		return false;

	if (t0.norm_square() < 1e-12 || t1.norm_square() < 1e-12)
		return false;

	double d = std::sqrt(p0.distance_square(p3));
	if (d < 1e-12)
		return false;

	double h = std::max(1e-6, d * 0.4);
	if (continuity == G1)
		return build_cubic_transition(p0, p3, t0, t1, h, out);

	if (continuity == G2)
		return build_quartic_fillet(p0, p3, t0, t1, d, h, out);

	return false;
}

bool NurbsCurveJoin::create_chamfer(const NurbsCurve& c1, const NurbsCurve& c2, double chamferLength, NurbsCurve& out)
{
	if (chamferLength <= 0 || c1.nb_points() < 2 || c2.nb_points() < 2)
		return false;

	Point3 p0 = c1.points().back();
	Point3 p3 = c2.points().front();
	Point3 t0, t1;
	if (!c1.tangent(1.0, t0) || !c2.tangent(0.0, t1))
		return false;

	if (t0.norm_square() < 1e-12 || t1.norm_square() < 1e-12)
		return false;

	double d = std::sqrt(p0.distance_square(p3));
	double usable = std::min(chamferLength, d * 0.49);
	if (usable <= 0)
		return false;

	Point3 ph0 = p0 + t0 * usable;
	Point3 ph1 = p3 - t1 * usable;
	return build_linear_connector(ph0, ph1, out);
}

bool NurbsCurveJoin::create_fillet(const NurbsCurve& c1, const NurbsCurve& c2, double radius, NurbsCurve& out)
{
	if (radius <= 0 || c1.nb_points() < 2 || c2.nb_points() < 2)
		return false;

	Point3 p0 = c1.points().back();
	Point3 p3 = c2.points().front();
	Point3 t0, t1;
	if (!c1.tangent(1.0, t0) || !c2.tangent(0.0, t1))
		return false;

	if (t0.norm_square() < 1e-12 || t1.norm_square() < 1e-12)
		return false;

	double d = std::sqrt(p0.distance_square(p3));
	if (d < 1e-12)
		return false;

	double h = std::min(radius, d * 0.45);
	return build_quartic_fillet(p0, p3, t0, t1, d, h, out);
}

