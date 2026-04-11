#include "NurbsBoolean.h"
#include "NurbsIntersection.h"

#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsTrimmedSurface.h"
#include "StepFile.h"
#include "Transform.h"
#include "OBJFile.h"
#include "NurbsUtil.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

void save_solid(const NurbsSolid& n, const string& filename)
{
	StepWriter sw;
	sw.open(filename + ".step");
	test_bool(sw.is_open(), string("step writer should open file: ") + filename);
	sw.write(n);
	sw.close();

	OBJWriter o;
	Mesh m;
	NurbsUtil::to_mesh(n, m);
	o.open(filename + ".obj");
	o.write(m);
}

void save_trimmed_surfaces(const std::vector<NurbsTrimmedSurface>& trimmed, const string& filename)
{
	StepWriter sw;
	sw.open(filename + ".step");
	test_bool(sw.is_open(), string("step writer should open file: ") + filename);
	for (const auto&  s :trimmed)
		sw.write(s);
	sw.close();

	OBJWriter o;
	Mesh m;
	NurbsUtil::to_mesh(trimmed, m, 26);
	o.open(filename + ".obj");
	o.write(m);
}

void check_trimmed_loop_invariants(const std::vector<NurbsTrimmedSurface>& trimmed, const string& tag)
{
	for (int i = 0; i < (int)trimmed.size(); ++i)
	{
		const auto& loops = trimmed[i].trim_loops();
		if (loops.empty())
			continue;

		int iOuter = 0;
		for (const auto& loop : loops)
		{
			test_bool(loop.points.size() >= 3, tag + ": trim loop should have at least 3 points");
			if (!loop.hole)
				iOuter++;
		}

		test_bool(iOuter > 0, tag + ": trimmed face should not contain hole-only loops");
	}
}

void test_nurbs_boolean_two_spheres_export()
{
	cout << endl << "test_nurbs_boolean_two_spheres_export" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_sphere(10., a);
	NurbsFactory::create_sphere(6., b);

	Translation t(Point3(25., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid u, i, d;

	const bool okU = nb.boolean_union(a, b, u);
	const bool okI = nb.boolean_intersection(a, b, i);
	const bool okD = nb.boolean_difference_bbox(a, b, d);

	test_bool(okU, "two-sphere union should succeed for disjoint solids");
	test_bool(okI, "two-sphere intersection should succeed for disjoint solids");
	test_bool(okD, "two-sphere difference should succeed for disjoint solids");
	test_bool(i.surfaces().empty(), "two-sphere disjoint intersection should be empty");

	save_solid(u, "test_nurbs_boolean_two_spheres_union");
	save_solid(d, "test_nurbs_boolean_two_spheres_diff");
}

void test_nurbs_boolean_disjoint()
{
	cout << endl << "test_nurbs_boolean_disjoint" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(4., 4., 4., a);
	NurbsFactory::create_box(2., 2., 2., b);
	Translation t(Point3(20., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid u, i, d;
	bool okU = nb.boolean_union(a, b, u);
	bool okI = nb.boolean_intersection(a, b, i);
	bool okD = nb.boolean_difference_bbox(a, b, d);

	test_bool(okU, "union should succeed on disjoint solids");
	test_bool(okI, "intersection should succeed on disjoint solids");
	test_bool(okD, "difference should succeed on disjoint solids");

	test_bool(u.surfaces().size() == a.surfaces().size() + b.surfaces().size(), "union disjoint surface count");
	test_bool(i.surfaces().empty(), "intersection disjoint should be empty");
	test_bool(d.surfaces().size() == a.surfaces().size(), "difference disjoint should keep A");
}

void test_nurbs_boolean_containment()
{
	cout << endl << "test_nurbs_boolean_containment" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(10., 10., 10., a);
	NurbsFactory::create_box(2., 2., 2., b);

	NurbsBoolean nb;
	NurbsSolid u, i, d;
	bool okU = nb.boolean_union(a, b, u);
	bool okI = nb.boolean_intersection(a, b, i);
	bool okD = nb.boolean_difference_bbox(a, b, d);

	test_bool(okU, "union should succeed for containment");
	test_bool(okI, "intersection should succeed for containment");
	test_bool(okD, "difference should succeed for containment");

	test_bool(u.surfaces().size() == a.surfaces().size(), "union containment should keep outer shell");
	test_bool(i.surfaces().size() == b.surfaces().size(), "intersection containment should keep inner shell");
	test_bool(d.surfaces().size() == a.surfaces().size() + b.surfaces().size(), "difference containment should add cavity shell");
}

void test_nurbs_boolean_box_cylinder_difference()
{
	cout << endl << "test_nurbs_boolean_box_cylinder_difference" << endl;

	NurbsSolid a, b, d;
	NurbsFactory::create_box(20., 20., 20., a);
	NurbsFactory::create_cylinder(6., 30., b);

	NurbsBoolean nb;
	bool ok = nb.compute_difference(a, b, d);
	test_bool(ok, "compute_difference should succeed for box minus cylinder");
	test_bool(!d.surfaces().empty(), "difference result should not be empty");
	test_bool(d.surfaces().size() > a.surfaces().size(), "difference should add cavity shell surfaces");
}

void test_nurbs_boolean_partial_overlap_is_reported()
{
	cout << endl << "test_nurbs_boolean_partial_overlap_is_reported" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(10., 10., 10., a);
	NurbsFactory::create_box(10., 10., 10., b);
	Translation t(Point3(5., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid u, i, d;

	test_bool(nb.boolean_union(a, b, u) == false, "union should report unsupported for partial overlap");
	test_bool(nb.boolean_intersection(a, b, i) == false, "intersection should report unsupported for partial overlap");
	test_bool(nb.boolean_difference_bbox(a, b, d) == false, "difference should report unsupported for partial overlap");
}

void test_nurbs_boolean_trimmed_pipeline_disjoint()
{
	cout << endl << "test_nurbs_boolean_trimmed_pipeline_disjoint" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(4., 4., 4., a);
	NurbsFactory::create_box(2., 2., 2., b);
	Translation t(Point3(20., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid resultSolid;
	NurbsIntersectionResult diagnostics;

	bool ok = nb.compute_union(a, b, resultSolid, &diagnostics);
	std::vector<NurbsTrimmedSurface> result;
	NurbsUtil::solid_to_trimmed_surfaces(resultSolid, result);
	test_bool(ok, "trimmed union should succeed for disjoint solids");
	test_bool(diagnostics.hasPartialOverlap == false, "diagnostics should not mark partial overlap");
	test_bool(result.size() == a.surfaces().size() + b.surfaces().size(), "trimmed disjoint union surface count");

	for (const auto& s : result)
		test_bool(!s.is_trimmed(), "disjoint case should not create trim loops");
}

void test_nurbs_boolean_trimmed_pipeline_partial_overlap_stub()
{
	cout << endl << "test_nurbs_boolean_trimmed_pipeline_partial_overlap_stub" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(10., 10., 10., a);
	NurbsFactory::create_box(10., 10., 10., b);
	Translation t(Point3(5., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid resultSolid;
	NurbsIntersectionResult diagnostics;

	bool ok = nb.compute_union(a, b, resultSolid, &diagnostics);
	std::vector<NurbsTrimmedSurface> result;
	if (ok)
		NurbsUtil::solid_to_trimmed_surfaces(resultSolid, result);
	test_bool(ok == false, "trimmed union should report unsupported overlap for now");
	test_bool(result.empty(), "trimmed result should be empty on unsupported overlap");
	test_bool(diagnostics.hasPartialOverlap, "diagnostics should flag partial overlap");
	test_bool(!diagnostics.curves.empty(), "diagnostics should provide stub curve");

	bool hasPolyline = false;
	bool checkedClosedWinding = false;
	for (const auto& c : diagnostics.curves)
	{
		if (c.samples.size() >= 2)
		{
			for (int i = 1; i < (int)c.samples.size(); ++i)
			{
				const Point3& p0 = c.samples[i - 1].point;
				const Point3& p1 = c.samples[i].point;
				test_bool(std::isfinite(p0.x()) && std::isfinite(p0.y()) && std::isfinite(p0.z()), "diagnostics point must be finite");
				test_bool(std::isfinite(p1.x()) && std::isfinite(p1.y()) && std::isfinite(p1.z()), "diagnostics point must be finite");
				test_bool((p1 - p0).norm_square() > 0., "diagnostics polyline should not contain duplicate adjacent points");
			}

			hasPolyline = true;

			if (c.closed && c.samples.size() >= 3)
			{
				double area2 = 0.;
				for (int i = 0; i < (int)c.samples.size(); ++i)
				{
					const auto& p0 = c.samples[i];
					const auto& p1 = c.samples[(i + 1) % c.samples.size()];
					area2 += p0.uA * p1.vA - p1.uA * p0.vA;
				}
				test_bool(area2 >= -1.e-12, "closed curve winding should be normalized (non-negative area in uvA)");
				checkedClosedWinding = true;
			}

			break;
		}
	}
	test_bool(hasPolyline, "diagnostics should contain at least one polyline curve");
	(void)checkedClosedWinding;
}

void test_nurbs_boolean_trimmed_intersection_partial_overlap_export()
{
	cout << endl << "test_nurbs_boolean_trimmed_intersection_partial_overlap_export" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_box(10., 10., 10., a);
	NurbsFactory::create_box(10., 10., 10., b);
	Translation t(Point3(5., 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid resultSolid;
	NurbsIntersectionResult diagnostics;

	bool ok = nb.compute_intersection(a, b, resultSolid, &diagnostics);
	std::vector<NurbsTrimmedSurface> result;
	if (ok)
		NurbsUtil::solid_to_trimmed_surfaces(resultSolid, result);
	test_bool(ok, "trimmed intersection should succeed for partial overlap approximation");
	test_bool(!result.empty(), "trimmed intersection should provide surfaces");
	test_bool(diagnostics.hasPartialOverlap, "diagnostics should flag partial overlap");
	test_bool(!diagnostics.curves.empty(), "diagnostics should contain intersection curves");

	save_trimmed_surfaces(result, "test_nurbs_boolean_trimmed_intersection_partial_overlap");

	OBJWriter ow;
	ow.open("test_nurbs_boolean_trimmed_intersection_partial_overlap_curves.obj");
	test_bool(ow.is_open(), "diagnostics curves obj should open");
	for (const auto& c : diagnostics.curves)
	{
		if (c.samples.size() < 2)
			continue;

		vector<Point3> polyline;
		polyline.reserve(c.samples.size());
		for (const auto& sample : c.samples)
			polyline.push_back(sample.point);

		ow.write(polyline);
	}
	ow.close();

	ifstream inCurves("test_nurbs_boolean_trimmed_intersection_partial_overlap_curves.obj");
	test_bool(inCurves.is_open(), "diagnostics curves obj should exist");
}

void test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap_export()
{
	cout << endl << "test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap_export" << endl;

	NurbsSolid torus, sphere;
	NurbsFactory::create_torus(8., 2., torus);
	NurbsFactory::create_sphere(8., sphere);

	NurbsBoolean nb;
	NurbsSolid classicIntersection;
	test_bool(nb.boolean_intersection(torus, sphere, classicIntersection) == false,
		"classic torus-sphere intersection should report unsupported overlap");

	std::vector<NurbsTrimmedSurface> result;
	NurbsIntersectionResult diagnostics;
	const bool ok = nb.compute_intersection(torus, sphere, classicIntersection, &diagnostics);
	
	if (ok)
		NurbsUtil::solid_to_trimmed_surfaces(classicIntersection, result);
	test_bool(ok, "trimmed torus-sphere intersection should succeed");
	test_bool(!result.empty(), "trimmed torus-sphere intersection should provide surfaces");
	test_bool(diagnostics.hasPartialOverlap, "diagnostics should flag partial overlap");
	test_bool(!diagnostics.curves.empty(), "diagnostics should contain torus-sphere intersection curves");

	save_trimmed_surfaces(result, "test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap");

	OBJWriter ow;
	ow.open("test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap_curves.obj");
	test_bool(ow.is_open(), "torus-sphere diagnostics curves obj should open");
	for (const auto& c : diagnostics.curves)
	{
		if (c.samples.size() < 2)
			continue;

		vector<Point3> polyline;
		polyline.reserve(c.samples.size());
		for (const auto& sample : c.samples)
			polyline.push_back(sample.point);

		ow.write(polyline);
	}
	ow.close();

	ifstream inCurves("test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap_curves.obj");
	test_bool(inCurves.is_open(), "torus-sphere diagnostics curves obj should exist");
}

void test_nurbs_boolean_exact_trimmed_torus_sphere_export()
{
	cout << endl << "test_nurbs_boolean_exact_trimmed_torus_sphere_export" << endl;

	NurbsSolid torus, sphere;
	NurbsFactory::create_torus(8., 2., torus);
	NurbsFactory::create_sphere(8., sphere);

	NurbsBoolean nb;

	NurbsSolid resUnionSolid, resInterSolid, resDiffSolid;
	NurbsIntersectionResult diagUnion, diagInter, diagDiff;

	const bool okU = nb.compute_union(torus, sphere, resUnionSolid, &diagUnion);
	std::vector<NurbsTrimmedSurface> resUnion;
	if (okU) NurbsUtil::solid_to_trimmed_surfaces(resUnionSolid, resUnion);
	const bool okI = nb.compute_intersection(torus, sphere, resInterSolid, &diagInter);
	std::vector<NurbsTrimmedSurface> resInter;
	if (okI) NurbsUtil::solid_to_trimmed_surfaces(resInterSolid, resInter);
	const bool okD = nb.compute_difference(torus, sphere, resDiffSolid, &diagDiff);
	std::vector<NurbsTrimmedSurface> resDiff;
	if (okD) NurbsUtil::solid_to_trimmed_surfaces(resDiffSolid, resDiff);

	test_bool(okU, "exact trimmed union should succeed on torus-sphere overlap");
	test_bool(okI, "exact trimmed intersection should succeed on torus-sphere overlap");
	test_bool(okD, "exact trimmed difference should succeed on torus-sphere overlap");

	test_bool(!resUnion.empty(), "exact trimmed union should produce surfaces");
	test_bool(!resInter.empty(), "exact trimmed intersection should produce surfaces");
	test_bool(!resDiff.empty(), "exact trimmed difference should produce surfaces");
	test_bool(resUnion.size() != resInter.size() || resDiff.size() != resInter.size(),
		"exact trimmed operations should not all collapse to the same surface count");

	test_bool(!diagUnion.curves.empty(), "exact trimmed union should collect intersection curves");
	test_bool(!diagInter.curves.empty(), "exact trimmed intersection should collect intersection curves");
	test_bool(!diagDiff.curves.empty(), "exact trimmed difference should collect intersection curves");

	check_trimmed_loop_invariants(resUnion, "exact torus-sphere union");
	check_trimmed_loop_invariants(resInter, "exact torus-sphere intersection");
	check_trimmed_loop_invariants(resDiff, "exact torus-sphere difference");

	save_trimmed_surfaces(resUnion, "test_nurbs_boolean_exact_trimmed_torus_sphere_union");
	save_trimmed_surfaces(resInter, "test_nurbs_boolean_exact_trimmed_torus_sphere_intersection");
	save_trimmed_surfaces(resDiff, "test_nurbs_boolean_exact_trimmed_torus_sphere_difference");
}

void test_nurbs_boolean_exact_trimmed_torus_torus_nested_loop_regression()
{
	cout << endl << "test_nurbs_boolean_exact_trimmed_torus_torus_nested_loop_regression" << endl;

	NurbsSolid a, b;
	NurbsFactory::create_torus(8., 2.5, a);
	NurbsFactory::create_torus(8., 2.5, b);

	RotationAngleAxis r(M_PI * 0.5, Point3(1., 0., 0.));
	b.apply_transform(r);
	Translation t(Point3(1.2, 0., 0.));
	b.apply_transform(t);

	NurbsBoolean nb;
	NurbsSolid resUnionSolid, resInterSolid, resDiffSolid;
	NurbsIntersectionResult diagUnion, diagInter, diagDiff;

	test_bool(nb.compute_union(a, b, resUnionSolid, &diagUnion), "exact torus-torus union should succeed");
	std::vector<NurbsTrimmedSurface> resUnion;
	if (nb.compute_union(a, b, resUnionSolid, &diagUnion)) NurbsUtil::solid_to_trimmed_surfaces(resUnionSolid, resUnion);
	test_bool(nb.compute_intersection(a, b, resInterSolid, &diagInter), "exact torus-torus intersection should succeed");
	std::vector<NurbsTrimmedSurface> resInter;
	if (nb.compute_intersection(a, b, resInterSolid, &diagInter)) NurbsUtil::solid_to_trimmed_surfaces(resInterSolid, resInter);
	test_bool(nb.compute_difference(a, b, resDiffSolid, &diagDiff), "exact torus-torus difference should succeed");
	std::vector<NurbsTrimmedSurface> resDiff;
	if (nb.compute_difference(a, b, resDiffSolid, &diagDiff)) NurbsUtil::solid_to_trimmed_surfaces(resDiffSolid, resDiff);

	test_bool(!diagUnion.curves.empty(), "exact torus-torus union should produce curves");
	test_bool(!diagInter.curves.empty(), "exact torus-torus intersection should produce curves");
	test_bool(!diagDiff.curves.empty(), "exact torus-torus difference should produce curves");

	check_trimmed_loop_invariants(resUnion, "exact torus-torus union");
	check_trimmed_loop_invariants(resInter, "exact torus-torus intersection");
	check_trimmed_loop_invariants(resDiff, "exact torus-torus difference");

	save_trimmed_surfaces(resInter, "test_nurbs_boolean_exact_trimmed_torus_torus_intersection");
}

void test_nurbs_boolean_random_cylinders_and_step_save()
{
	cout << endl << "test_nurbs_boolean_random_cylinders_and_step_save" << endl;

	srand(1337);
	NurbsBoolean nb;

	for (int iCase = 0; iCase < 5; ++iCase)
	{
		const double dRadiusA = 1. + 4. * (rand() / (double)RAND_MAX);
		const double dRadiusB = 1. + 4. * (rand() / (double)RAND_MAX);
		const double dHeightA = 2. + 8. * (rand() / (double)RAND_MAX);
		const double dHeightB = 2. + 8. * (rand() / (double)RAND_MAX);

		NurbsSolid a, b;
		NurbsFactory::create_cylinder(dRadiusA, dHeightA, a);
		NurbsFactory::create_cylinder(dRadiusB, dHeightB, b);

		const double dMargin = 2. + 4. * (rand() / (double)RAND_MAX);
		Translation t(Point3(dRadiusA + dRadiusB + dMargin, 0., 0.));
		b.apply_transform(t);

		NurbsSolid u, inter, diff;
		const bool okU = nb.boolean_union(a, b, u);
		const bool okI = nb.boolean_intersection(a, b, inter);
		const bool okD = nb.boolean_difference_bbox(a, b, diff);

		test_bool(okU, "random cylinder union should succeed for disjoint solids");
		test_bool(okI, "random cylinder intersection should succeed for disjoint solids");
		test_bool(okD, "random cylinder difference should succeed for disjoint solids");
		test_bool(inter.surfaces().empty(), "random disjoint cylinder intersection should be empty");

		const string sBase = string("test_random_cylinder_bool_") + to_string(iCase);
		save_solid(a, sBase + "_a.step");
		save_solid(b, sBase + "_b.step");
		save_solid(u, sBase + "_union.step");
		save_solid(diff, sBase + "_diff.step");
	}
}

int main()
{
	test_nurbs_boolean_two_spheres_export();
	test_nurbs_boolean_disjoint();
	//test_nurbs_boolean_containment();
	test_nurbs_boolean_box_cylinder_difference();
	test_nurbs_boolean_partial_overlap_is_reported();
	test_nurbs_boolean_trimmed_pipeline_disjoint();
	//test_nurbs_boolean_trimmed_pipeline_partial_overlap_stub();
	test_nurbs_boolean_trimmed_intersection_partial_overlap_export();
	//test_nurbs_boolean_trimmed_intersection_torus_sphere_overlap_export();
	//test_nurbs_boolean_mesh_fallback_torus_sphere_export();
	//test_nurbs_boolean_exact_trimmed_torus_sphere_export();
	test_nurbs_boolean_exact_trimmed_torus_torus_nested_loop_regression();
	test_nurbs_boolean_random_cylinders_and_step_save();

	cout << "Test Finished.";
	return 0;
}
