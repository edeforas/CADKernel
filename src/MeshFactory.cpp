#include "MeshFactory.h"

#include "MeshTessellate.h"
#include "MeshUtil.h"

#include <cmath>
using namespace std;

#include "Transform.h"

namespace MeshFactory
{
	void create_tetrahedron(double dSize, Mesh& m)
	{
		m.clear();

		m.add_vertex(dSize, dSize, dSize);
		m.add_vertex(-dSize, -dSize, dSize);
		m.add_vertex(dSize, -dSize, -dSize);
		m.add_vertex(-dSize, dSize, -dSize);

		m.add_triangle(0, 1, 2);
		m.add_triangle(1, 3, 2);
		m.add_triangle(0, 3, 1);
		m.add_triangle(0, 2, 3);
	}

	void create_box(double xSize, double ySize, double zSize, Mesh& m)
	{
		double x = xSize / 2.;
		double y = ySize / 2.;
		double z = zSize / 2.;

		m.clear();

		m.add_vertex(x, y, z);
		m.add_vertex(-x, y, z);
		m.add_vertex(-x, -y, z);
		m.add_vertex(x, -y, z);
		m.add_vertex(x, y, -z);
		m.add_vertex(-x, y, -z);
		m.add_vertex(-x, -y, -z);
		m.add_vertex(x, -y, -z);

		m.add_quad(0, 1, 2, 3);
		m.add_quad(7, 6, 5, 4);
		m.add_quad(4, 5, 1, 0);
		m.add_quad(3, 2, 6, 7);
		m.add_quad(5, 6, 2, 1);
		m.add_quad(0, 3, 7, 4);
	}
	///////////////////////////////////////////////////////////////////////////
	void create_octahedron(double dSize, Mesh& m)
	{
		m.clear();

		m.add_vertex(dSize, 0, 0);
		m.add_vertex(0., dSize, 0);
		m.add_vertex(0., 0., dSize);
		m.add_vertex(-dSize, 0, 0);
		m.add_vertex(0, -dSize, 0.);
		m.add_vertex(0., 0, -dSize);

		m.add_triangle(0, 1, 2);
		m.add_triangle(1, 0, 5);
		m.add_triangle(2, 1, 3);
		m.add_triangle(0, 2, 4);
		m.add_triangle(4, 3, 5);
		m.add_triangle(3, 4, 2);
		m.add_triangle(5, 3, 1);
		m.add_triangle(4, 5, 0);
	}

	void create_dodecahedron(double dSize, Mesh& m)
	{
		m.clear();
		// coordinate from: https://en.wikipedia.org/wiki/Regular_dodecahedron

		double u = dSize;
		double t1 = (std::sqrt(5.) + 1.) / 2.;
		double t2 = 1. / t1;
		t1 *= dSize;
		t2 *= dSize;

		m.add_vertex(u, u, u);
		m.add_vertex(u, u, -u);
		m.add_vertex(u, -u, u);
		m.add_vertex(u, -u, -u);
		m.add_vertex(-u, u, u);
		m.add_vertex(-u, u, -u);
		m.add_vertex(-u, -u, u);
		m.add_vertex(-u, -u, -u);

		m.add_vertex(0, t1, t2);
		m.add_vertex(0, t1, -t2);
		m.add_vertex(0, -t1, t2);
		m.add_vertex(0, -t1, -t2);

		m.add_vertex(t2, 0, t1);
		m.add_vertex(t2, 0, -t1);
		m.add_vertex(-t2, 0, t1);
		m.add_vertex(-t2, 0, -t1);

		m.add_vertex(t1, t2, 0);
		m.add_vertex(t1, -t2, 0);
		m.add_vertex(-t1, t2, 0);
		m.add_vertex(-t1, -t2, 0);

		m.add_pentagon(0, 12, 2, 17, 16);
		m.add_pentagon(19, 18, 5, 15, 7);

		m.add_pentagon(12, 0, 8, 4, 14);
		m.add_pentagon(13, 3, 11, 7, 15);

		m.add_pentagon(2, 12, 14, 6, 10);
		m.add_pentagon(9, 1, 13, 15, 5);

		m.add_pentagon(18, 19, 6, 14, 4);
		m.add_pentagon(3, 13, 1, 16, 17);

		m.add_pentagon(1, 9, 8, 0, 16);
		m.add_pentagon(19, 7, 11, 10, 6);

		m.add_pentagon(18, 4, 8, 9, 5);
		m.add_pentagon(2, 10, 11, 3, 17);
	}

	void create_icosahedron(double dSize, Mesh& m)
	{
		m.clear();

		// from "Beyond the Third Dimension" by Banchoff, p.168.
		double u = dSize;
		double t = (std::sqrt(5.) - 1.) / 2. * dSize;

		m.add_vertex(u, 0., t); //0
		m.add_vertex(u, 0., -t);
		m.add_vertex(-u, 0., t);
		m.add_vertex(-u, 0., -t);

		m.add_vertex(0., t, u); //4
		m.add_vertex(0., -t, u);
		m.add_vertex(0., t, -u);
		m.add_vertex(0., -t, -u);

		m.add_vertex(t, u, 0.); //8
		m.add_vertex(-t, u, 0.);
		m.add_vertex(t, -u, 0.);
		m.add_vertex(-t, -u, 0.);

		m.add_triangle(0, 1, 8);
		m.add_triangle(0, 10, 1);
		m.add_triangle(2, 3, 11);
		m.add_triangle(2, 9, 3);

		m.add_triangle(4, 5, 0);
		m.add_triangle(4, 2, 5);
		m.add_triangle(6, 7, 3);
		m.add_triangle(6, 1, 7);

		m.add_triangle(8, 9, 4);
		m.add_triangle(8, 6, 9);
		m.add_triangle(10, 11, 7);
		m.add_triangle(10, 5, 11);

		m.add_triangle(0, 8, 4);
		m.add_triangle(3, 7, 11);

		m.add_triangle(10, 0, 5);
		m.add_triangle(9, 6, 3);

		m.add_triangle(1, 6, 8);
		m.add_triangle(2, 11, 5);

		m.add_triangle(2, 4, 9);
		m.add_triangle(1, 10, 7);
	}

	void create_cylinder(double dHeight, double dRadius, int iNbSegments, Mesh& m)
	{
		m.clear();

		double PI = 3.14159265358979323846264338; // TODO
		double z = dHeight / 2.;

		int iNbSegmentLongitude = 2 + 4 * iNbSegments;
		// create vertices
		for (int i = 0; i < iNbSegmentLongitude; i++)
		{
			double x = std::cos((double)i / iNbSegmentLongitude * 2. * PI) * dRadius;
			double y = std::sin((double)i / iNbSegmentLongitude * 2. * PI) * dRadius;

			m.add_vertex(x, y, z);
			m.add_vertex(x, y, -z);
		}
		m.add_vertex(0, 0, z);
		m.add_vertex(0, 0, -z);

		// associates triangles, edge 
		for (int i = 0; i < iNbSegmentLongitude - 1; i++)
		{
			m.add_quad(2 * i, 2 * i + 1, 2 * i + 3, 2 * i + 2);
		}
		m.add_quad(2 * iNbSegmentLongitude - 2, 2 * iNbSegmentLongitude - 1, 1, 0);

		// associates triangles, up face 
		for (int i = 0; i < iNbSegmentLongitude; i++)
		{
			m.add_triangle(iNbSegmentLongitude * 2, 2 * i, 2 * i + 2);
		}
		m.add_triangle(iNbSegmentLongitude * 2, 2 * iNbSegmentLongitude - 2, 0);

		// associates triangles, down face 
		for (int i = 0; i < iNbSegmentLongitude; i++)
		{
			m.add_triangle(iNbSegmentLongitude * 2 + 1, 2 * i + 3, 2 * i + 1);
		}
		m.add_triangle(iNbSegmentLongitude * 2 + 1, 1, 2 * iNbSegmentLongitude - 1);
	}

	void create_sphere_geodesic(double dRadius, int iNbSegments, Mesh& m)
	{
		m.clear();

		Mesh mi;
		create_icosahedron(1., mi);

		MeshTessellate mt;
		mt.compute(mi, iNbSegments, m);

		for (int i = 0; i < m.nb_vertices(); i++)
		{
			Point3 v;
			m.get_vertex(i, v);
			m.set_vertex(i, v.normalized() * dRadius);
		}

		double dEdgesize = dRadius / 0.951; // as in https://en.wikipedia.org/wiki/Regular_icosahedron

		MeshUtil::merge_vertices(m,dEdgesize / iNbSegments / 20.); //todo better formula
	}

	void create_sphere_uv(double dRadius, int iNbSegments, Mesh& m)
	{
		m.clear();

		double PI = 3.14159265358979323846264338; // TODO

		int iNbSegmentLongitude = 2 + 4 * iNbSegments; //at least an octahedron
		int iNbSegmentLatitude = 2 + 4 * iNbSegments / 2; //at least an octahedron

		// create vertices
		for (int longitude = 0; longitude < iNbSegmentLongitude; longitude++)
		{
			double upi = (double)longitude / iNbSegmentLongitude * 2. * PI;

			for (int latitude = 1; latitude < iNbSegmentLatitude; latitude++)
			{
				double vpi = (double)latitude / iNbSegmentLatitude * PI;

				// compute parametric point from: https://en.wikipedia.org/wiki/Sphere
				Point3 p;
				p._x = dRadius * std::sin(vpi) * std::cos(upi);
				p._y = dRadius * std::sin(vpi) * std::sin(upi);
				p._z = dRadius * std::cos(vpi);

				m.add_vertex(p);
			}
		}
		m.add_vertex(0, 0, dRadius); // pole north
		m.add_vertex(0, 0, -dRadius); //pole south

		// associates quads to vertices
		int iNbVertexLatitude = iNbSegmentLatitude - 1; //poles removed
		//int iNbVertex = m.nb_vertices()-2; // poles removed
		int iPoleNorth = m.nb_vertices() - 2;
		int iPoleSouth = m.nb_vertices() - 1;

		for (int longitude = 0; longitude < iNbSegmentLongitude; longitude++)
		{
			// add first trig
			int iLongitude1 = longitude;
			int iLongitude2 = (longitude + 1) % iNbSegmentLongitude;

			// make north pole
			m.add_triangle(iPoleNorth, iLongitude1 * iNbVertexLatitude, iLongitude2 * iNbVertexLatitude);

			for (int latitude = 1; latitude < iNbVertexLatitude; latitude++)
			{
				m.add_quad(
					iLongitude1 * iNbVertexLatitude + latitude - 1,
					iLongitude1 * iNbVertexLatitude + latitude,
					iLongitude2 * iNbVertexLatitude + latitude,
					iLongitude2 * iNbVertexLatitude + latitude - 1
				);
			}

			// make south pole
			m.add_triangle(iLongitude1 * iNbVertexLatitude + iNbVertexLatitude - 1, iPoleSouth, iLongitude2 * iNbVertexLatitude + iNbVertexLatitude - 1);
		}
	}

	void create_torus(double dMajorRadius, double dMinorRadius, int iNbSegments, Mesh& m)
	{
		m.clear();

		if (!m.empty())
			return;

		double PI = 3.14159265358979323846264338; // TODO

		int iNbSegmentLongitude = 2 + 4 * iNbSegments;
		int iNbSegmentLatitude = 2 + 4 + 4 * iNbSegments / 2;

		// create vertices
		for (int longitude = 0; longitude < iNbSegmentLongitude; longitude++)
		{
			double upi = (double)longitude / iNbSegmentLongitude * 2. * PI;

			for (int latitude = 0; latitude < iNbSegmentLatitude; latitude++)
			{
				double vpi = (double)latitude / iNbSegmentLatitude * 2. * PI;

				// compute parametric point from: https://en.wikipedia.org/wiki/Torus
				Point3 p;
				p._x = (dMajorRadius + dMinorRadius * std::cos(vpi)) * std::cos(upi);
				p._y = (dMajorRadius + dMinorRadius * std::cos(vpi)) * std::sin(upi);
				p._z = dMinorRadius * std::sin(vpi);

				m.add_vertex(p);
			}
		}

		// associates quads
		for (int longitude = 0; longitude < iNbSegmentLongitude; longitude++)
		{
			for (int latitude = 0; latitude < iNbSegmentLatitude; latitude++)
			{
				int longitude2 = (longitude + 1) % iNbSegmentLongitude;
				int latitude2 = (latitude + 1) % iNbSegmentLatitude;
				m.add_quad(
					longitude2 * iNbSegmentLatitude + latitude,
					longitude2 * iNbSegmentLatitude + latitude2,
					longitude * iNbSegmentLatitude + latitude2,
					longitude * iNbSegmentLatitude + latitude);
			}
		}
	}

	void create_revolve(const std::vector<Point3>& profile, int iNbSegments, Mesh& m)
	{
		create_revolve(profile, iNbSegments, false, m);
	}

	void create_revolve(const std::vector<Point3>& profile, int iNbSegments, bool bCapEnds, Mesh& m)
	{
		m.clear();

		if ((int)profile.size() < 2)
			return;

		if (iNbSegments < 3)
			iNbSegments = 3;

		double PI = 3.14159265358979323846264338;
		const int iNbProfilePoints = (int)profile.size();

		for (int segment = 0; segment < iNbSegments; ++segment)
		{
			double angle = (double)segment / iNbSegments * 2. * PI;
			double c = std::cos(angle);
			double s = std::sin(angle);

			for (int i = 0; i < iNbProfilePoints; ++i)
			{
				const Point3& p = profile[i];
				m.add_vertex(
					p.x() * c - p.y() * s,
					p.x() * s + p.y() * c,
					p.z());
			}
		}

		for (int segment = 0; segment < iNbSegments; ++segment)
		{
			const int segmentNext = (segment + 1) % iNbSegments;

			for (int i = 0; i < iNbProfilePoints - 1; ++i)
			{
				m.add_quad(
					segmentNext * iNbProfilePoints + i,
					segmentNext * iNbProfilePoints + i + 1,
					segment * iNbProfilePoints + i + 1,
					segment * iNbProfilePoints + i);
			}
		}

		if (!bCapEnds)
			return;

		const int iFirst = 0;
		const int iLast = iNbProfilePoints - 1;
		const Point3& pFirst = profile[iFirst];
		const Point3& pLast = profile[iLast];
		const double eps = 1.e-12;

		const double rFirst = std::sqrt(pFirst.x() * pFirst.x() + pFirst.y() * pFirst.y());
		if (rFirst > eps)
		{
			const int iCenterFirst = m.add_vertex(0., 0., pFirst.z());
			for (int segment = 0; segment < iNbSegments; ++segment)
			{
				const int segmentNext = (segment + 1) % iNbSegments;
				const int iV1 = segment * iNbProfilePoints + iFirst;
				const int iV2 = segmentNext * iNbProfilePoints + iFirst;
				m.add_triangle(iCenterFirst, iV2, iV1);
			}
		}

		const double rLast = std::sqrt(pLast.x() * pLast.x() + pLast.y() * pLast.y());
		if (rLast > eps)
		{
			const int iCenterLast = m.add_vertex(0., 0., pLast.z());
			for (int segment = 0; segment < iNbSegments; ++segment)
			{
				const int segmentNext = (segment + 1) % iNbSegments;
				const int iV1 = segment * iNbProfilePoints + iLast;
				const int iV2 = segmentNext * iNbProfilePoints + iLast;
				m.add_triangle(iCenterLast, iV1, iV2);
			}
		}
	}
}








