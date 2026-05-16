#include "NurbsUtil.h"
#include "NurbsBasis.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsFactory.h"
#include "BezierSurface.h"
#include "Mesh.h"

#include <cmath>
#include <algorithm>

namespace
{
	bool point_in_uv_loop(const std::vector<NurbsUvPoint>& loop, double u, double v)
	{
		if (loop.size() < 3)
			return false;

		bool inside = false;
		for (int i = 0, j = (int)loop.size() - 1; i < (int)loop.size(); j = i++)
		{
			const double xi = loop[i].u;
			const double yi = loop[i].v;
			const double xj = loop[j].u;
			const double yj = loop[j].v;

			const bool intersect = ((yi > v) != (yj > v)) &&
				(u < (xj - xi) * (v - yi) / ((yj - yi) + 1.e-20) + xi);

			if (intersect)
				inside = !inside;
		}

		return inside;
	}

	bool point_inside_trims(const NurbsTrimmedSurface& ts, double u, double v)
	{
		if (!ts.is_trimmed())
			return true;

		bool insideOuter = false;
		bool insideHole = false;

		for (const auto& loop : ts.trim_loops())
		{
			if (!point_in_uv_loop(loop.points, u, v))
				continue;

			if (loop.hole)
				insideHole = true;
			else
				insideOuter = true;
		}

		return insideOuter && !insideHole;
	}

	static BezierSurface extract_bezier_patch(const NurbsSurface& n, size_t uSpan, size_t vSpan,
	                                         const NurbsUtil::KnotAnalysis& uKnots, const NurbsUtil::KnotAnalysis& vKnots);
	static int find_span_start_index(const std::vector<double>& knots, double knotValue, int degree);

	// Helper function to extract a single Bezier patch from a NURBS surface span
	static BezierSurface extract_bezier_patch(const NurbsSurface& n, size_t uSpan, size_t vSpan,
	                                         const NurbsUtil::KnotAnalysis& uKnots, const NurbsUtil::KnotAnalysis& vKnots)
	{
		BezierSurface patch;

		// Create a copy of the surface to work with
		NurbsSurface workSurface = n;

		// Insert knots to achieve Bezier multiplicity for this span
		// For Bezier form, we need multiplicity = degree + 1 at each knot
		int bezierMultiplicityU = n.degree_u() + 1;
		int bezierMultiplicityV = n.degree_v() + 1;

		// Insert knots at the span boundaries to achieve Bezier multiplicity
		double uStart = uKnots.unique_knots[uSpan];
		double uEnd = uKnots.unique_knots[uSpan + 1];
		double vStart = vKnots.unique_knots[vSpan];
		double vEnd = vKnots.unique_knots[vSpan + 1];

		// Insert knots to achieve Bezier multiplicity at boundaries
		for (int i = uKnots.multiplicities[uSpan]; i < bezierMultiplicityU; ++i) {
			workSurface.insert_knot_u(uStart);
		}
		for (int i = uKnots.multiplicities[uSpan + 1]; i < bezierMultiplicityU; ++i) {
			workSurface.insert_knot_u(uEnd);
		}

		for (int i = vKnots.multiplicities[vSpan]; i < bezierMultiplicityV; ++i) {
			workSurface.insert_knot_v(vStart);
		}
		for (int i = vKnots.multiplicities[vSpan + 1]; i < bezierMultiplicityV; ++i) {
			workSurface.insert_knot_v(vEnd);
		}

		// Find the control points for this span
		// After knot insertion, the control points for the span are consecutive
		int uStartIndex = find_span_start_index(workSurface.knots_u(), uStart, n.degree_u());
		int vStartIndex = find_span_start_index(workSurface.knots_v(), vStart, n.degree_v());

		// Extract (degreeU+1) x (degreeV+1) control points
		std::vector<Point3> bezierPoints;
		for (int i = 0; i <= n.degree_u(); ++i) {
			for (int j = 0; j <= n.degree_v(); ++j) {
				int pointIndex = (uStartIndex + i) * workSurface.nb_points_v() + (vStartIndex + j);
				if (pointIndex < (int)workSurface.points().size()) {
					bezierPoints.push_back(workSurface.points()[pointIndex]);
				}
			}
		}

		// Create the Bezier patch
		if ((int)bezierPoints.size() == (n.degree_u() + 1) * (n.degree_v() + 1)) {
			patch.set_degree(n.degree_u(), n.degree_v());
			patch.set_control_points(bezierPoints, n.degree_u() + 1, n.degree_v() + 1);
		}

		return patch;
	}

	// Helper function to find the starting index of a knot span
	static int find_span_start_index(const std::vector<double>& knots, double knotValue, int degree)
	{
		for (size_t i = degree; i < knots.size() - degree; ++i) {
			if (std::abs(knots[i] - knotValue) < 1e-10) {
				return i - degree;
			}
		}
		return 0; // Default to start
	}
}

///////////////////////////////////////////////////////////////////////////
void NurbsUtil::create_curve_from_points(const std::vector<Point3>& points, int iDegree, NurbsCurve& n)
{
	n.clear();

	if (points.empty())
		return;

	int iMaxDegree = (int)points.size() - 1;
	if (iMaxDegree < 0)
		iMaxDegree = 0;

	int iFinalDegree = iDegree;
	if (iFinalDegree < 0)
		iFinalDegree = 0;
	if (iFinalDegree > iMaxDegree)
		iFinalDegree = iMaxDegree;

	n.set_degree(iFinalDegree);
	n.set_points(points);
	n.set_uniform();
	n.set_equals_weights();
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::create_surface_from_z(const std::vector<double>& z, int iSizeX, int iSizeY, int iDegree, NurbsSurface& n)
{
	n.clear();

	std::vector < Point3> points;

	int idx = 0;
	for (int x = 0; x < iSizeX; x++)
		for (int y = 0; y < iSizeY; y++)
		{
			Point3 p;
			p.x() = x;
			p.y() = y;
			p.z() = z[idx];
			points.push_back(p);

			idx++;
		}

	n.set_degree(iDegree, iDegree);
	n.set_points(points, iSizeX, iSizeY);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::create_solid_from_mesh(const Mesh& m, NurbsSolid& n)
{
	n.clear();
	for (int i = 0; i < m.nb_triangles(); i++)
	{
		Triangle3 t;
		m.get_triangle(i, t);

		NurbsSurface nsurf;
		NurbsFactory::create_triangle(t.p1(), t.p2(), t.p3(), nsurf);

		n.add_surface(nsurf);
	}
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::to_controlpoints_mesh(const NurbsSurface& n, Mesh& m) //show the ctrl points mesh lattice
{
	const std::vector<Point3>& points = n.points();
	int iNbPointsU = n.nb_points_u();
	int iNbPointsV = n.nb_points_v();

	m.clear();
	for (int u = 0; u < iNbPointsU - 1; u++)
		for (int v = 0; v < iNbPointsV - 1; v++)
		{
			Point3 p1, p2, p3, p4;

			p1 = points[v * iNbPointsU + u];
			p2 = points[v * iNbPointsU + (u + 1)];
			p3 = points[(v + 1) * iNbPointsU + (u + 1)];
			p4 = points[(v + 1) * iNbPointsU + u];

			m.add_quad(p1, p2, p3, p4);
		}
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::to_mesh(const NurbsSurface& n, Mesh& m, int iNbSegments, bool bClearMesh)
{
	if (bClearMesh)
		m.clear();

	bool bIsClosedU = n.is_closed_u();
	bool bIsClosedV = n.is_closed_v();

	int iNbPointsStart = m.nb_vertices();
	int iNbPointsU = iNbSegments * n.nb_points_u();
	int iNbPointsV = iNbSegments * n.nb_points_v();

	if ((iNbPointsU == 0) || (iNbPointsV == 0))
		return;

	// add vertices, for now we keep all vertices, even if not used because closed
	Point3 p;
	for (int v = 0; v <= iNbPointsV; v++)
		for (int u = 0; u <= iNbPointsU; u++)
		{
			double du1 = (double)u / iNbPointsU;
			double dv1 = (double)v / iNbPointsV;

			n.evaluate(du1, dv1, p);
			m.add_vertex(p);
		}

	// add quad linked to vertices
	for (int v = 0; v < iNbPointsV; v++)
		for (int u = 0; u < iNbPointsU; u++)
		{
			int iEndU = u + 1;
			int iEndV = v + 1;

			if (bIsClosedU && (u == iNbPointsU - 1))
				iEndU = 0;

			if (bIsClosedV && (v == iNbPointsV - 1))
				iEndV = 0;

			m.add_quad(
				iNbPointsStart + u + (iNbPointsU + 1) * v,
				iNbPointsStart + iEndU + (iNbPointsU + 1) * v,
				iNbPointsStart + iEndU + (iNbPointsU + 1) * iEndV,
				iNbPointsStart + u + (iNbPointsU + 1) * iEndV
			);
		}
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::to_mesh(const NurbsSolid& ns, Mesh& m, int iNbSegments)
{
	m.clear();
	for (const auto& f : ns.surfaces())
	{
		to_mesh(f, m, iNbSegments, false);
	}
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::to_mesh(const NurbsTrimmedSurface& ts, Mesh& m, int iNbSegments, bool bClearMesh)
{
	if (iNbSegments < 2)
		iNbSegments = 2;

	if (bClearMesh)
		m.clear();

	const NurbsSurface& s = ts;
	const int stride = iNbSegments + 1;
	const int iStartVertex = m.nb_vertices();

	for (int iv = 0; iv <= iNbSegments; ++iv)
	{
		double v = (double)iv / (double)iNbSegments;
		for (int iu = 0; iu <= iNbSegments; ++iu)
		{
			double u = (double)iu / (double)iNbSegments;
			Point3 p;
			s.evaluate(u, v, p);
			m.add_vertex(p);
		}
	}

	for (int iv = 0; iv < iNbSegments; ++iv)
	{
		for (int iu = 0; iu < iNbSegments; ++iu)
		{
			double uc = ((double)iu + 0.5) / (double)iNbSegments;
			double vc = ((double)iv + 0.5) / (double)iNbSegments;
			if (!point_inside_trims(ts, uc, vc))
				continue;

			const int i00 = iStartVertex + iu + stride * iv;
			const int i10 = iStartVertex + (iu + 1) + stride * iv;
			const int i11 = iStartVertex + (iu + 1) + stride * (iv + 1);
			const int i01 = iStartVertex + iu + stride * (iv + 1);

			m.add_triangle(i00, i10, i11);
			m.add_triangle(i00, i11, i01);
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void NurbsUtil::to_mesh(const std::vector<NurbsTrimmedSurface>& trimmedSurfaces, Mesh& m, int iNbSegments)
{
	m.clear();
	for (int i = 0; i < (int)trimmedSurfaces.size(); i++)
		to_mesh(trimmedSurfaces[i], m, iNbSegments, false);
}
///////////////////////////////////////////////////////////////////////////
double NurbsUtil::sanitize_weight(double value, double kEpsilonWeight)
{
	if (!std::isfinite(value))
		return 1.;

	value = std::fabs(value);
	if (value < kEpsilonWeight)
		return kEpsilonWeight;

	return value;
}
///////////////////////////////////////////////////////////////////////////
std::vector<double> NurbsUtil::build_safe_weights(const std::vector<double>& weights, int expectedSize)
{
	std::vector<double> safeWeights;
	safeWeights.reserve(expectedSize);

	if ((int)weights.size() != expectedSize)
	{
		safeWeights.assign(expectedSize, 1.);
		return safeWeights;
	}

	for (int i = 0; i < expectedSize; ++i)
		safeWeights.push_back(NurbsUtil::sanitize_weight(weights[i]));

	return safeWeights;
}
///////////////////////////////////////////////////////////////////////////
std::vector<double> NurbsUtil::build_segmented_quadratic_knots(int nbSegments)
{
	return NurbsBasis::build_segmented_quadratic_knots(nbSegments);
}
///////////////////////////////////////////////////////////////////////////
std::vector<double> NurbsUtil::build_uniform_knots(int degree, int nbCtrlPoints)
{
	return NurbsBasis::build_uniform_knots(degree, nbCtrlPoints);
}


NurbsUtil::KnotAnalysis NurbsUtil::analyze_knots(const std::vector<double>& knots, int expectedDegree, int expectedCtrlPoints)
{
	KnotAnalysis result;

	std::vector<double> work = knots;
	const int expectedCount = expectedCtrlPoints + expectedDegree + 1;
	if (expectedCtrlPoints <= 0)
		return result;

	if ((int)work.size() != expectedCount)
		work = NurbsUtil::build_uniform_knots(expectedDegree, expectedCtrlPoints);

	if (work.empty())
		return result;

	for (auto& k : work)
	{
		if (!std::isfinite(k))
			k = 0;
	}

	std::sort(work.begin(), work.end());

	const double eps = 1.e-9;
	double current = work[0];
	int count = 1;
	for (int i = 1; i < (int)work.size(); ++i)
	{
		if (std::fabs(work[i] - current) <= eps)
		{
			count++;
		}
		else
		{
			result.unique_knots.push_back(current);
			result.multiplicities.push_back(count);
			current = work[i];
			count = 1;
		}
	}

	result.unique_knots.push_back(current);
	result.multiplicities.push_back(count);

	return result;
}

void NurbsUtil::to_bezier_patches(const NurbsSurface& n, std::vector<BezierSurface>& patches)
{
	patches.clear();

	// Analyze knot vectors to find spans
	KnotAnalysis uKnots = analyze_knots(n.knots_u(), n.degree_u(), n.nb_points_u());
	KnotAnalysis vKnots = analyze_knots(n.knots_v(), n.degree_v(), n.nb_points_v());

	if (uKnots.unique_knots.size() < 2 || vKnots.unique_knots.size() < 2) {
		return; // Not enough spans
	}

	// For each span [u_i, u_{i+1}] x [v_j, v_{j+1}], create a Bezier patch
	for (size_t ui = 0; ui < uKnots.unique_knots.size() - 1; ++ui) {
		for (size_t vi = 0; vi < vKnots.unique_knots.size() - 1; ++vi) {
			BezierSurface patch = extract_bezier_patch(n, ui, vi, uKnots, vKnots);
			if (patch.is_valid()) {
				patches.push_back(patch);
			}
		}
	}
}

void NurbsUtil::solid_to_trimmed_surfaces(const NurbsSolid& src, std::vector<NurbsTrimmedSurface>& dst)
{
	dst.clear();
	dst.reserve(src.surfaces().size());

	for (const auto& s : src.surfaces())
	{
		NurbsTrimmedSurface ts(s);
		dst.push_back(ts);
	}
}