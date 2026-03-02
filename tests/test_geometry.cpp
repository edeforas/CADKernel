#include "Geometry.h"

#include <string>
#include <iostream>
using namespace std;

void test(bool b, string s)
{
	if (b == false)
	{
		cout << "condition: [" << s << "] is not realized!" << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_plane_segment_intersection()
{
	Point3 p1(1., 0., 0.);
	Point3 p2(0., 1., 0.);
	Point3 p3(0., 0., 1.);
	Plane3 plane1(p1, p2, p3);

	Point3 p4(0., 0., 0.); // cube diagonal
	Point3 p5(1., 1., 1.);
	Point3 intersection;

	bool bSucceded = plane1.intersect_with(Segment3(p4, p5), intersection);
	test(bSucceded, "Segment must intersect");
	test((intersection - Point3(1., 1., 1.) / 3.).norm() < 1.e-10, "Intersection must be at (1 / 3, 1 / 3, 1 / 3)");
}

///////////////////////////////////////////////////////////////////////////
void test_plane_point_projection()
{
	Point3 p1(1., 0., 0.);
	Point3 p2(0., 1., 0.);
	Point3 p3(0., 0., 1.);
	Plane3 plane1(p1, p2, p3);

	Point3 p5(1., 1., 1.);
	Point3 projected;

	plane1.project_point(p5, projected);

	test((projected - Point3(1., 1., 1.) / 3.).norm() < 1.e-10, "Projected must be at (1/3,1/3,1/3)");
}
///////////////////////////////////////////////////////////////////////////
void test_point_in_triangle()
{
	Point3 p1(1., 0., 0.);
	Point3 p2(0., 1., 0.);
	Point3 p3(0., 0., 1.);
	Triangle3 t3(p1, p2, p3);

	Point3 pA(0.33, 0.33, 0.33);
	test(t3.contains(pA), "Point A is inside");

	Point3 pB(1., 1., -1.); //out of plane
	test(t3.contains(pB) == false, "Point B is outside");
}
///////////////////////////////////////////////////////////////////////////
void test_triangle_triangle_intersection()
{
	// triangles nearly touching at the (0,0,0) point
	{
		Point3 pA1(0., 1., 0.);
		Point3 pA2(0., -1., 0.);
		Point3 pA3(-1., 0., 0.);
		Triangle3 tA(pA1, pA2, pA3);

		Point3 pB1(0.1, 0., 1.);
		Point3 pB2(0.1, 0., -1.);
		Point3 pB3(1., 0., 0.);
		Triangle3 tB(pB1, pB2, pB3);

		test(tA.intersect_with(tB) == false, "TriangleA and triangleB must not intersect");
	}
	{
		Point3 pA1(0., 1., 0.);
		Point3 pA2(0., -1., 0.);
		Point3 pA3(-1., 0., 0.);
		Triangle3 tA(pA1, pA2, pA3);

		Point3 pB1(-0.1, 0., 1.);
		Point3 pB2(-0.1, 0., -1.);
		Point3 pB3(1., 0., 0.);
		Triangle3 tB(pB1, pB2, pB3);

		test(tA.intersect_with(tB), "TriangleA and triangleB must intersect");
	}

	// coplanar overlap (z=0): should intersect
	{
		Triangle3 tA(Point3(0., 0., 0.), Point3(2., 0., 0.), Point3(0., 2., 0.));
		Triangle3 tB(Point3(1., 0., 0.), Point3(3., 0., 0.), Point3(1., 2., 0.));

		test(tA.intersect_with(tB), "Coplanar overlapping triangles must intersect");
	}

	// coplanar disjoint (z=0): should not intersect
	{
		Triangle3 tA(Point3(0., 0., 0.), Point3(1., 0., 0.), Point3(0., 1., 0.));
		Triangle3 tB(Point3(2., 2., 0.), Point3(3., 2., 0.), Point3(2., 3., 0.));

		test(tA.intersect_with(tB) == false, "Coplanar disjoint triangles must not intersect");
	}

	// coplanar touching at one vertex: should intersect
	{
		Triangle3 tA(Point3(0., 0., 0.), Point3(1., 0., 0.), Point3(0., 1., 0.));
		Triangle3 tB(Point3(1., 0., 0.), Point3(2., 0., 0.), Point3(1., 1., 0.));

		test(tA.intersect_with(tB), "Coplanar triangles touching at a vertex must intersect");
	}

	// coplanar sharing one full edge: should intersect
	{
		Triangle3 tA(Point3(0., 0., 0.), Point3(1., 0., 0.), Point3(0., 1., 0.));
		Triangle3 tB(Point3(1., 0., 0.), Point3(0., 0., 0.), Point3(1., -1., 0.));

		test(tA.intersect_with(tB), "Coplanar triangles sharing an edge must intersect");
	}
}

///////////////////////////////////////////////////////////////////////////
void test_segment_intersection()
{
	Point3 p1(1., 0., 0.);
	Point3 p2(0., 1., 0.);
	Point3 p3(-1., 0., 0.);
	Point3 p4(0., -1., 0.);

	Point3 pInter(10., 20., 30.);

	bool bIntersect = Segment3(p1, p3).intersect(Segment3(p2, p4), pInter);
	test(bIntersect, "Segments must intersect");
	test((pInter - Point3(0., 0., 0.)).norm() < 1.e-10, "Intersection must be at origin");

	bIntersect = Segment3(p1, p2).intersect(Segment3(p3, p4), pInter);
	test(bIntersect == false, "Segments must not intersect");
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	cout << "Tests start..." << endl;

	test_plane_point_projection();
	test_plane_segment_intersection();
	test_point_in_triangle();
	test_triangle_triangle_intersection();
	test_segment_intersection();

	cout << "Tests finished." << endl;
	return 0;
}