#include "MeshBoolean.h"
#include "MeshFactory.h"
#include "Transform.h"
#include "OBJFile.h"
#include "STLFile.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

void test(bool a, const string& sMessage = "")
{
	if (a == false)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

int nb_linked_triangles(const Mesh& m)
{
	int iCount = 0;
	for (int i = 0; i < m.nb_triangles(); ++i)
	{
		if (!m.is_triangle_unlinked(i))
			iCount++;
	}
	return iCount;
}

int nb_vertices_on_both_spheres(const Mesh& m, const Point3& centerA, double radiusA, const Point3& centerB, double radiusB, double tol)
{
	int iCount = 0;
	for (int i = 0; i < m.nb_vertices(); ++i)
	{
		Point3 p;
		m.get_vertex(i, p);

		double dA = std::fabs((p - centerA).norm() - radiusA);
		double dB = std::fabs((p - centerB).norm() - radiusB);
		if ((dA <= tol) && (dB <= tol))
			iCount++;
	}

	return iCount;
}

void test_nested_boxes_split()
{
	cout << endl << "test_nested_boxes_split" << endl;

	Mesh boxOuter;
	Mesh boxInner;
	MeshFactory::create_box(4., 4., 4., boxOuter);
	MeshFactory::create_box(2., 2., 2., boxInner);

	MeshBoolean mb;
	Mesh Aonly, Bonly, AInB, BInA;
	mb.compute_split(boxOuter, boxInner, Aonly, Bonly, AInB, BInA);

	test(nb_linked_triangles(Aonly) == 12, "outer box should remain outside");
	test(nb_linked_triangles(Bonly) == 0, "inner box should be fully inside outer box");
	test(nb_linked_triangles(AInB) == 0, "outer box should not be inside inner box");
	test(nb_linked_triangles(BInA) == 12, "inner box triangles should be classified inside outer box");

	Mesh intersection;
	mb.compute_intersection(boxOuter, boxInner, intersection);
	test(nb_linked_triangles(intersection) == 12, "nested boxes intersection should equal inner box shell");
}

void test_partial_overlap_boxes_split()
{
	cout << endl << "test_partial_overlap_boxes_split" << endl;

	Mesh boxA;
	Mesh boxB;
	MeshFactory::create_box(2., 2., 2., boxA);
	MeshFactory::create_box(2., 2., 2., boxB);

	Translation tr(Point3(1., 0., 0.));
	boxB.apply_transform(tr);

	MeshBoolean mb;
	Mesh Aonly, Bonly, AInB, BInA;
	mb.compute_split(boxA, boxB, Aonly, Bonly, AInB, BInA);

	test(nb_linked_triangles(Aonly) > 0, "A must keep triangles outside B");
	test(nb_linked_triangles(Bonly) > 0, "B must keep triangles outside A");
	test(nb_linked_triangles(AInB) > 0, "A must have triangles classified inside B");
	test(nb_linked_triangles(BInA) > 0, "B must have triangles classified inside A");

	Mesh intersection;
	mb.compute_intersection(boxA, boxB, intersection);
	test(nb_linked_triangles(intersection) > 0, "partial overlap intersection must not be empty");
}

void test_union_and_difference_nested_boxes()
{
	cout << endl << "test_union_and_difference_nested_boxes" << endl;

	Mesh boxOuter;
	Mesh boxInner;
	MeshFactory::create_box(4., 4., 4., boxOuter);
	MeshFactory::create_box(2., 2., 2., boxInner);

	MeshBoolean mb;

	Mesh meshUnion;
	mb.compute_union(boxOuter, boxInner, meshUnion);
	test(nb_linked_triangles(meshUnion) == 12, "nested boxes union should equal outer box shell");

	Mesh meshDifference;
	mb.compute_difference(boxOuter, boxInner, meshDifference);
	test(nb_linked_triangles(meshDifference) == 24, "outer-inner difference should contain outer shell and inner cavity shell");
}

void test_difference_two_spheres()
{
	cout << endl << "test_difference_two_spheres" << endl;

	Mesh sphereA;
	Mesh sphereB;
	MeshFactory::create_sphere_uv(1., 8, sphereA);
	MeshFactory::create_sphere_uv(1., 8, sphereB);

	Translation tr(Point3(0.9, 0., 0.));
	sphereB.apply_transform(tr);

	MeshBoolean mb;
	Mesh difference;
	mb.compute_difference(sphereA, sphereB, difference);

	int nbA = nb_linked_triangles(sphereA);
	int nbDiff = nb_linked_triangles(difference);

	test(nbDiff > 0, "sphere subtraction should not be empty");
	test(nbDiff > nbA, "sphere subtraction should add inner cut shell triangles");

	Mesh intersection;
	mb.compute_intersection(sphereA, sphereB, intersection);
	test(nb_linked_triangles(intersection) > 0, "overlapping spheres should have non-empty intersection");

	const Point3 centerA(0., 0., 0.);
	const Point3 centerB(0.9, 0., 0.);
	const double radius = 1.;
	const double tol = 5.e-2;

	int iSeamDiff = nb_vertices_on_both_spheres(difference, centerA, radius, centerB, radius, tol);
	int iSeamInter = nb_vertices_on_both_spheres(intersection, centerA, radius, centerB, radius, tol);
	test(iSeamDiff >= 8, "sphere subtraction should contain seam vertices lying on both spheres");
	test(iSeamInter >= 8, "sphere intersection should contain seam vertices lying on both spheres");

	Mesh meshUnion;
	mb.compute_union(sphereA, sphereB, meshUnion);

	bool bSaved = OBJFile::save("test_meshboolean_spheres_AminusB.obj", difference);
	test(bSaved, "failed to save sphere subtraction OBJ");

	bSaved = OBJFile::save("test_meshboolean_spheres_intersection.obj", intersection);
	test(bSaved, "failed to save sphere intersection OBJ");

	bSaved = OBJFile::save("test_meshboolean_spheres_union.obj", meshUnion);
	test(bSaved, "failed to save sphere union OBJ");
}

void test_cube_intersection_stl()
{
	cout << endl << "test_cube_intersection_stl" << endl;

	Mesh cubeA;
	Mesh cubeB;
	MeshFactory::create_box(2., 2., 2., cubeA);
	MeshFactory::create_box(2., 2., 2., cubeB);

	Translation tr(Point3(0.8, 0.5, 0.3));
	cubeB.apply_transform(tr);

	MeshBoolean mb;
	Mesh intersection;
	mb.compute_intersection(cubeA, cubeB, intersection);
	Mesh meshUnion;
	mb.compute_union(cubeA, cubeB, meshUnion);
	Mesh meshDifference;
	mb.compute_difference(cubeA, cubeB, meshDifference);

	test(nb_linked_triangles(intersection) > 0, "cube intersection should not be empty");
	test(nb_linked_triangles(meshUnion) > 0, "cube union should not be empty");
	test(nb_linked_triangles(meshDifference) > 0, "cube difference should not be empty");

	bool bSaved = STLFile::save("test_meshboolean_cubes_intersection.stl", intersection);
	test(bSaved, "failed to save cube intersection STL");

	bSaved = STLFile::save("test_meshboolean_cubes_union.stl", meshUnion);
	test(bSaved, "failed to save cube union STL");

	bSaved = STLFile::save("test_meshboolean_cubes_AminusB.stl", meshDifference);
	test(bSaved, "failed to save cube difference STL");
}

int main()
{
	test_nested_boxes_split();
	test_partial_overlap_boxes_split();
	test_union_and_difference_nested_boxes();
	test_difference_two_spheres();
	test_cube_intersection_stl();

	cout << "Test Finished.";
	return 0;
}
