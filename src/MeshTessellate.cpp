#include "MeshTessellate.h"

///////////////////////////////////////////////////////////////////////////
MeshTessellate::MeshTessellate()
{
}

MeshTessellate::~MeshTessellate()
{
}

void MeshTessellate::compute(const Mesh& meshIn, int iNbSegments, Mesh& meshOut)
{
	meshOut.clear();
	for (int iTriangle = 0; iTriangle < meshIn.nb_triangles(); iTriangle++)
	{
		if (meshIn.is_triangle_unlinked(iTriangle))
			continue;

		Point3 a, b, c;
		int iA, iB, iC;

		meshIn.get_triangle_vertices(iTriangle, iA, iB, iC);
		meshIn.get_vertex(iA, a);
		meshIn.get_vertex(iB, b);
		meshIn.get_vertex(iC, c);
		//meshOut.add_triangle(a, b, c);

		//direct tessellating by adding small stripes of triangles
		for (int stripe = 0; stripe < iNbSegments; stripe++)
		{
			//add first triangle
			Point3 p1x = interpolate(a, b, c, stripe, 0, iNbSegments);
			Point3 p2x = interpolate(a, b, c, stripe + 1, 0, iNbSegments);
			Point3 p3x = interpolate(a, b, c, stripe, 1, iNbSegments);
			meshOut.add_triangle(p1x, p2x, p3x);

			//finish stripe
			for (int j = 0; j < stripe; j++)
			{
				Point3 p1 = interpolate(a, b, c, stripe - j, j, iNbSegments);
				Point3 p2 = interpolate(a, b, c, stripe - j, j + 1, iNbSegments);
				Point3 p3 = interpolate(a, b, c, stripe - j - 1, j + 1, iNbSegments);
				Point3 p4 = interpolate(a, b, c, stripe - j - 1, j + 2, iNbSegments);
				meshOut.add_triangle(p1, p2, p3);
				meshOut.add_triangle(p3, p2, p4);
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////
Point3 MeshTessellate::interpolate(const Point3& a, const Point3& b, const Point3& c, int uStep, int vStep, int maxStep)
{
	double tu = (double)uStep / (double)maxStep;
	double tv = (double)vStep / (double)maxStep;

	return a + (b - a) * tu + (c - a) * tv;
}