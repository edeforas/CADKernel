#include "NurbsFillet.h"

#include "NurbsCurve.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
	Point3 estimate_tangent(const std::vector<Point3>& points, int i)
	{
		const int n = (int)points.size();
		if (n <= 1)
			return Point3(1., 0., 0.);

		Point3 t;
		if (i <= 0)
			t = points[1] - points[0];
		else if (i >= n - 1)
			t = points[n - 1] - points[n - 2];
		else
			t = points[i + 1] - points[i - 1];

		if (t.norm_square() < 1.e-20)
			return Point3(1., 0., 0.);
		return t.normalized();
	}



	bool extract_edge_curve(const NurbsSurface& s, NurbsFillet::SurfaceEdge edge, NurbsCurve& c)
	{
		c.clear();

		const int nU = s.nb_points_u();
		const int nV = s.nb_points_v();
		if (nU <= 0 || nV <= 0)
			return false;

		const std::vector<Point3>& pts = s.points();
		const std::vector<double>& w = s.weights();
		if ((int)pts.size() != nU * nV || (int)w.size() != nU * nV)
			return false;

		std::vector<Point3> outPts;
		std::vector<double> outW;
		int degree = 1;
		std::vector<double> knots;

		if ((edge == NurbsFillet::EdgeVMin) || (edge == NurbsFillet::EdgeVMax))
		{
			const int v = (edge == NurbsFillet::EdgeVMin) ? 0 : (nV - 1);
			outPts.reserve(nU);
			outW.reserve(nU);
			for (int u = 0; u < nU; ++u)
			{
				const int idx = v * nU + u;
				outPts.push_back(pts[idx]);
				outW.push_back(w[idx]);
			}
			degree = s.degree_u();
			knots = s.knots_u();
		}
		else
		{
			const int u = (edge == NurbsFillet::EdgeUMin) ? 0 : (nU - 1);
			outPts.reserve(nV);
			outW.reserve(nV);
			for (int v = 0; v < nV; ++v)
			{
				const int idx = v * nU + u;
				outPts.push_back(pts[idx]);
				outW.push_back(w[idx]);
			}
			degree = s.degree_v();
			knots = s.knots_v();
		}

		c.set_degree(std::max(1, degree));
		c.set_points(outPts);
		c.set_weights(outW);
		c.set_knots(knots);
		return true;
	}

	double curve_pair_error(const NurbsCurve& c1, const NurbsCurve& c2)
	{
		if (c1.nb_points() <= 0 || c1.nb_points() != c2.nb_points())
			return 1.e300;

		const std::vector<Point3>& p1 = c1.points();
		const std::vector<Point3>& p2 = c2.points();

		double errForward = 0.;
		double errReverse = 0.;
		for (int i = 0; i < (int)p1.size(); ++i)
		{
			errForward += p1[i].distance_square(p2[i]);
			errReverse += p1[i].distance_square(p2[(int)p2.size() - 1 - i]);
		}

		return std::min(errForward, errReverse) / (double)p1.size();
	}
}

NurbsFillet::NurbsFillet()
{
}

NurbsFillet::~NurbsFillet()
{
}

bool NurbsFillet::create_chamfer(const NurbsCurve& c1, const NurbsCurve& c2, double dOffset1, double dOffset2, NurbsSurface& out) const
{
	out.clear();

	if (!NurbsCurveUtil::curves_are_compatible(c1, c2))
		return false;
	if (!std::isfinite(dOffset1) || !std::isfinite(dOffset2))
		return false;
	if (dOffset1 < 0. || dOffset2 < 0.)
		return false;

	const std::vector<Point3>& p1 = c1.points();
	const std::vector<Point3>& p2 = c2.points();
	const std::vector<double>& w1 = c1.weights();
	const std::vector<double>& w2 = c2.weights();

	std::vector<Point3> pts;
	std::vector<double> w;
	pts.reserve(2 * p1.size());
	w.reserve(2 * p1.size());

	for (int i = 0; i < (int)p1.size(); ++i)
	{
		const Point3 d = p2[i] - p1[i];
		const double len = d.norm();
		if (len < 1.e-12)
		{
			pts.push_back(p1[i]);
			w.push_back(w1[i]);
			pts.push_back(p2[i]);
			w.push_back(w2[i]);
			continue;
		}

		double t1 = dOffset1 / len;
		double t2 = dOffset2 / len;
		t1 = std::max(0., std::min(0.49, t1));
		t2 = std::max(0., std::min(0.49, t2));
		if (t1 + t2 > 0.98)
		{
			const double s = 0.98 / (t1 + t2);
			t1 *= s;
			t2 *= s;
		}

		const Point3 q1 = p1[i] + d * t1;
		const Point3 q2 = p2[i] - d * t2;

		pts.push_back(q1);
		w.push_back(std::max(1.e-12, std::fabs(w1[i])));
		pts.push_back(q2);
		w.push_back(std::max(1.e-12, std::fabs(w2[i])));
	}

	out.set_degree(c1.degree(), 1);
	out.set_points(pts, c1.nb_points(), 2);
	out.set_knots_u(c1.knots());
	out.set_knots_v({ 0., 0., 1., 1. });
	out.set_weights(w);
	out.set_closed_u(c1.is_closed() && c2.is_closed());
	out.set_closed_v(false);
	return true;
}

bool NurbsFillet::create_fillet(const NurbsCurve& c1, const NurbsCurve& c2, double dRadius, NurbsSurface& out) const
{
	out.clear();

	if (!NurbsCurveUtil::curves_are_compatible(c1, c2))
		return false;
	if (!std::isfinite(dRadius) || dRadius < 0.)
		return false;

	const std::vector<Point3>& p1 = c1.points();
	const std::vector<Point3>& p2 = c2.points();
	const std::vector<double>& w1 = c1.weights();
	const std::vector<double>& w2 = c2.weights();

	std::vector<Point3> pts;
	std::vector<double> w;
	pts.reserve(3 * p1.size());
	w.reserve(3 * p1.size());

	for (int i = 0; i < (int)p1.size(); ++i)
	{
		const Point3 d = p2[i] - p1[i];
		const double len = d.norm();
		if (len < 1.e-12)
		{
			pts.push_back(p1[i]);
			w.push_back(std::max(1.e-12, std::fabs(w1[i])));
			pts.push_back(p1[i]);
			w.push_back(std::max(1.e-12, std::fabs(w1[i])));
			pts.push_back(p2[i]);
			w.push_back(std::max(1.e-12, std::fabs(w2[i])));
			continue;
		}

		double t = dRadius / len;
		t = std::max(0., std::min(0.45, t));

		const Point3 q1 = p1[i] + d * t;
		const Point3 q2 = p2[i] - d * t;

		const Point3 tan1 = estimate_tangent(p1, i);
		const Point3 tan2 = estimate_tangent(p2, i);
		Point3 tan = tan1 + tan2;
		if (tan.norm_square() < 1.e-20)
			tan = tan1;
		if (tan.norm_square() < 1.e-20)
			tan = Point3(1., 0., 0.);

		Point3 n = tan.cross_product(d);
		if (n.norm_square() < 1.e-20)
			n = Point3(0., 0., 1.).cross_product(d);
		if (n.norm_square() < 1.e-20)
			n = Point3(0., 1., 0.);
		n.normalize();

		const Point3 mid = (q1 + q2) * 0.5;
		const Point3 qm = mid + n * dRadius;

		const double ww1 = std::max(1.e-12, std::fabs(w1[i]));
		const double ww2 = std::max(1.e-12, std::fabs(w2[i]));
		const double wwm = std::max(1.e-12, std::sqrt(ww1 * ww2) * std::sqrt(0.5));

		pts.push_back(q1);
		w.push_back(ww1);
		pts.push_back(qm);
		w.push_back(wwm);
		pts.push_back(q2);
		w.push_back(ww2);
	}

	out.set_degree(c1.degree(), 2);
	out.set_points(pts, c1.nb_points(), 3);
	out.set_knots_u(c1.knots());
	out.set_knots_v({ 0., 0., 0., 1., 1., 1. });
	out.set_weights(w);
	out.set_closed_u(c1.is_closed() && c2.is_closed());
	out.set_closed_v(false);
	return true;
}

bool NurbsFillet::create_chamfer_between_surfaces(const NurbsSurface& s1, SurfaceEdge e1, const NurbsSurface& s2, SurfaceEdge e2, double dOffset1, double dOffset2, NurbsSurface& out) const
{
	out.clear();

	NurbsCurve c1, c2;
	if (!extract_edge_curve(s1, e1, c1))
		return false;
	if (!extract_edge_curve(s2, e2, c2))
		return false;
	if (!NurbsCurveUtil::align_curve_orientation(c1, c2))
		return false;

	return create_chamfer(c1, c2, dOffset1, dOffset2, out);
}

bool NurbsFillet::create_fillet_between_surfaces(const NurbsSurface& s1, SurfaceEdge e1, const NurbsSurface& s2, SurfaceEdge e2, double dRadius, NurbsSurface& out) const
{
	out.clear();

	NurbsCurve c1, c2;
	if (!extract_edge_curve(s1, e1, c1))
		return false;
	if (!extract_edge_curve(s2, e2, c2))
		return false;
	if (!NurbsCurveUtil::align_curve_orientation(c1, c2))
		return false;

	return create_fillet(c1, c2, dRadius, out);
}

bool NurbsFillet::find_shared_edge_pair(const NurbsSurface& s1, const NurbsSurface& s2, SurfaceEdge& e1, SurfaceEdge& e2, double dTol) const
{
	const SurfaceEdge edges[4] = { EdgeUMin, EdgeUMax, EdgeVMin, EdgeVMax };

	double bestErr = 1.e300;
	SurfaceEdge bestE1 = EdgeUMin;
	SurfaceEdge bestE2 = EdgeUMin;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			NurbsCurve c1, c2;
			if (!extract_edge_curve(s1, edges[i], c1))
				continue;
			if (!extract_edge_curve(s2, edges[j], c2))
				continue;

			if (c1.nb_points() != c2.nb_points())
				continue;

			const double err = curve_pair_error(c1, c2);
			if (err < bestErr)
			{
				bestErr = err;
				bestE1 = edges[i];
				bestE2 = edges[j];
			}
		}
	}

	if (!std::isfinite(bestErr))
		return false;
	if (bestErr > dTol * dTol)
		return false;

	e1 = bestE1;
	e2 = bestE2;
	return true;
}

bool NurbsFillet::create_chamfer_on_solid(const NurbsSolid& solid, int iSurfaceA, SurfaceEdge eA, int iSurfaceB, SurfaceEdge eB, double dOffset1, double dOffset2, NurbsSurface& out) const
{
	out.clear();

	const std::vector<NurbsSurface>& surfaces = solid.surfaces();
	if (iSurfaceA < 0 || iSurfaceB < 0)
		return false;
	if (iSurfaceA >= (int)surfaces.size() || iSurfaceB >= (int)surfaces.size())
		return false;

	return create_chamfer_between_surfaces(surfaces[iSurfaceA], eA, surfaces[iSurfaceB], eB, dOffset1, dOffset2, out);
}

bool NurbsFillet::create_fillet_on_solid(const NurbsSolid& solid, int iSurfaceA, SurfaceEdge eA, int iSurfaceB, SurfaceEdge eB, double dRadius, NurbsSurface& out) const
{
	out.clear();

	const std::vector<NurbsSurface>& surfaces = solid.surfaces();
	if (iSurfaceA < 0 || iSurfaceB < 0)
		return false;
	if (iSurfaceA >= (int)surfaces.size() || iSurfaceB >= (int)surfaces.size())
		return false;

	return create_fillet_between_surfaces(surfaces[iSurfaceA], eA, surfaces[iSurfaceB], eB, dRadius, out);
}

bool NurbsFillet::create_chamfer_on_solid_auto(const NurbsSolid& solid, int iSurfaceA, int iSurfaceB, double dOffset1, double dOffset2, NurbsSurface& out) const
{
	out.clear();

	const std::vector<NurbsSurface>& surfaces = solid.surfaces();
	if (iSurfaceA < 0 || iSurfaceB < 0)
		return false;
	if (iSurfaceA >= (int)surfaces.size() || iSurfaceB >= (int)surfaces.size())
		return false;

	SurfaceEdge eA = EdgeUMin;
	SurfaceEdge eB = EdgeUMin;
	if (!find_shared_edge_pair(surfaces[iSurfaceA], surfaces[iSurfaceB], eA, eB))
		return false;

	return create_chamfer_between_surfaces(surfaces[iSurfaceA], eA, surfaces[iSurfaceB], eB, dOffset1, dOffset2, out);
}

bool NurbsFillet::create_fillet_on_solid_auto(const NurbsSolid& solid, int iSurfaceA, int iSurfaceB, double dRadius, NurbsSurface& out) const
{
	out.clear();

	const std::vector<NurbsSurface>& surfaces = solid.surfaces();
	if (iSurfaceA < 0 || iSurfaceB < 0)
		return false;
	if (iSurfaceA >= (int)surfaces.size() || iSurfaceB >= (int)surfaces.size())
		return false;

	SurfaceEdge eA = EdgeUMin;
	SurfaceEdge eB = EdgeUMin;
	if (!find_shared_edge_pair(surfaces[iSurfaceA], surfaces[iSurfaceB], eA, eB))
		return false;

	return create_fillet_between_surfaces(surfaces[iSurfaceA], eA, surfaces[iSurfaceB], eB, dRadius, out);
}

bool NurbsFillet::create_fillet_or_chamfer_on_solid(const NurbsSolid& solid, int iSurfaceA, int iSurfaceB, bool isFillet, double size, NurbsSurface& out) const
{
	if (isFillet)
	{
		return create_fillet_on_solid_auto(solid, iSurfaceA, iSurfaceB, size, out);
	}
	else
	{
		return create_chamfer_on_solid_auto(solid, iSurfaceA, iSurfaceB, size, size, out);
	}
}
