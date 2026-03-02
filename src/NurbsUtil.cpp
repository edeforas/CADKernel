#include "NurbsUtil.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsFactory.h"
#include "Mesh.h"

#include <cmath>

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
void NurbsUtil::create_from_z(const std::vector<double>& z, int iSizeX, int iSizeY, int iDegree, NurbsSurface& n)
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
void NurbsUtil::create_from_mesh(const Mesh& m, NurbsSolid& n)
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
double NurbsUtil::sanitize_weight(double value,double kEpsilonWeight )
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
	std::vector<double> knots;
	if (nbSegments <= 0)
		return knots;

	knots.reserve(2 * nbSegments + 4);
	knots.push_back(0.);
	knots.push_back(0.);
	knots.push_back(0.);

	for (int i = 1; i < nbSegments; ++i)
	{
		const double t = (double)i / (double)nbSegments;
		knots.push_back(t);
		knots.push_back(t);
	}

	knots.push_back(1.);
	knots.push_back(1.);
	knots.push_back(1.);
	return knots;
}
///////////////////////////////////////////////////////////////////////////
