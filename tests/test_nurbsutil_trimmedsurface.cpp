#include "NurbsFactory.h"
#include "NurbsSurface.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsUtil.h"
#include "Mesh.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

void test_nurbsutil_trimmedsurface_outer_clip()
{
	cout << endl << "test_nurbsutil_trimmedsurface_outer_clip" << endl;

	NurbsSurface base;
	NurbsFactory::create_quad(
		Point3(0., 0., 0.),
		Point3(1., 0., 0.),
		Point3(1., 1., 0.),
		Point3(0., 1., 0.),
		base);

	NurbsTrimmedSurface ts(base);
	ts.add_outer_loop({
		{0.2, 0.2},
		{0.8, 0.2},
		{0.8, 0.8},
		{0.2, 0.8}
	});

	Mesh fullMesh;
	NurbsUtil::to_mesh(base, fullMesh, 20);

	Mesh trimmedMesh;
	NurbsUtil::to_mesh(ts, trimmedMesh, 20);

	test_bool(fullMesh.nb_triangles() > 0, "full mesh should contain triangles");
	test_bool(trimmedMesh.nb_triangles() > 0, "trimmed mesh should contain triangles");
	test_bool(trimmedMesh.nb_triangles() < fullMesh.nb_triangles(), "outer clipping should reduce triangle count");
}

void test_nurbsutil_trimmedsurface_hole_clip()
{
	cout << endl << "test_nurbsutil_trimmedsurface_hole_clip" << endl;

	NurbsSurface base;
	NurbsFactory::create_quad(
		Point3(0., 0., 0.),
		Point3(1., 0., 0.),
		Point3(1., 1., 0.),
		Point3(0., 1., 0.),
		base);

	NurbsTrimmedSurface ts(base);
	ts.add_outer_loop({
		{0.0, 0.0},
		{1.0, 0.0},
		{1.0, 1.0},
		{0.0, 1.0}
	});
	ts.add_inner_loop({
		{0.35, 0.35},
		{0.65, 0.35},
		{0.65, 0.65},
		{0.35, 0.65}
	});

	Mesh fullMesh;
	NurbsUtil::to_mesh(base, fullMesh, 24);

	Mesh holeMesh;
	NurbsUtil::to_mesh(ts, holeMesh, 24);

	test_bool(holeMesh.nb_triangles() > 0, "hole-trimmed mesh should contain triangles");
	test_bool(holeMesh.nb_triangles() < fullMesh.nb_triangles(), "hole clipping should reduce triangle count");
}

void test_nurbsutil_trimmedsurface_vector_overload()
{
	cout << endl << "test_nurbsutil_trimmedsurface_vector_overload" << endl;

	NurbsSurface left, right;
	NurbsFactory::create_quad(
		Point3(0., 0., 0.),
		Point3(1., 0., 0.),
		Point3(1., 1., 0.),
		Point3(0., 1., 0.),
		left);
	NurbsFactory::create_quad(
		Point3(1., 0., 0.),
		Point3(2., 0., 0.),
		Point3(2., 1., 0.),
		Point3(1., 1., 0.),
		right);

	NurbsTrimmedSurface tsLeft(left), tsRight(right);
	tsLeft.add_outer_loop({ {0.1, 0.1}, {0.9, 0.1}, {0.9, 0.9}, {0.1, 0.9} });
	tsRight.add_outer_loop({ {0.1, 0.1}, {0.9, 0.1}, {0.9, 0.9}, {0.1, 0.9} });

	Mesh mLeft, mRight, mBoth;
	NurbsUtil::to_mesh(tsLeft, mLeft, 18);
	NurbsUtil::to_mesh(tsRight, mRight, 18);
	NurbsUtil::to_mesh(std::vector<NurbsTrimmedSurface>{ tsLeft, tsRight }, mBoth, 18);

	test_bool(mLeft.nb_triangles() > 0 && mRight.nb_triangles() > 0, "single trimmed meshes should contain triangles");
	test_bool(mBoth.nb_triangles() == mLeft.nb_triangles() + mRight.nb_triangles(), "vector overload should aggregate meshes");
}

double signed_area_uv(const std::vector<NurbsUvPoint>& loop)
{
	double area2 = 0.;
	for (int i = 0; i < (int)loop.size(); ++i)
	{
		const auto& p0 = loop[i];
		const auto& p1 = loop[(i + 1) % loop.size()];
		area2 += p0.u * p1.v - p1.u * p0.v;
	}
	return area2 * 0.5;
}

void test_nurbsutil_trimmedsurface_loop_autofix()
{
	cout << endl << "test_nurbsutil_trimmedsurface_loop_autofix" << endl;

	NurbsSurface base;
	NurbsFactory::create_quad(
		Point3(0., 0., 0.),
		Point3(1., 0., 0.),
		Point3(1., 1., 0.),
		Point3(0., 1., 0.),
		base);

	NurbsTrimmedSurface ts(base);

	// outer loop with duplicates, explicit closure and clockwise order
	ts.add_outer_loop({
		{0.2, 0.2}, {0.2, 0.2}, {0.2, 0.8}, {0.8, 0.8}, {0.8, 0.2}, {0.2, 0.2}
	});

	// inner loop with duplicates, explicit closure and counter-clockwise order (should be flipped)
	ts.add_inner_loop({
		{0.35, 0.35}, {0.65, 0.35}, {0.65, 0.65}, {0.35, 0.65}, {0.35, 0.65}, {0.35, 0.35}
	});

	test_bool(ts.trim_loops().size() == 2, "two valid loops should be kept after normalization");

	const auto& outer = ts.trim_loops()[0];
	const auto& inner = ts.trim_loops()[1];

	test_bool(!outer.hole, "first loop should be outer");
	test_bool(inner.hole, "second loop should be hole");
	test_bool(outer.points.size() >= 3, "outer loop should have valid polygon size");
	test_bool(inner.points.size() >= 3, "inner loop should have valid polygon size");

	test_bool((outer.points.front().u != outer.points.back().u) || (outer.points.front().v != outer.points.back().v), "outer loop should be stored without duplicated closing point");
	test_bool((inner.points.front().u != inner.points.back().u) || (inner.points.front().v != inner.points.back().v), "inner loop should be stored without duplicated closing point");

	test_bool(signed_area_uv(outer.points) > 0., "outer winding should be normalized CCW");
	test_bool(signed_area_uv(inner.points) < 0., "inner winding should be normalized CW");

	Mesh m;
	NurbsUtil::to_mesh(ts, m, 22);
	test_bool(m.nb_triangles() > 0, "autofixed trimmed loop should still generate mesh");
}

int main()
{
	test_nurbsutil_trimmedsurface_outer_clip();
	test_nurbsutil_trimmedsurface_hole_clip();
	test_nurbsutil_trimmedsurface_vector_overload();
	test_nurbsutil_trimmedsurface_loop_autofix();

	cout << "Test Finished.";
	return 0;
}
