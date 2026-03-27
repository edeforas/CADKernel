#include "NurbsExtrude.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "NurbsFactory.h"

#include "OBJFile.h"
#include "StepFile.h"

#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

void test_near(double a, double ref, double epsilon=1.e-10,const string& sMessage="")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsextrude_cylinder()
{
	cout << endl << "test_nurbsextrude_cylinder" << endl;

	//create profile curve
	NurbsCurve nc;
	NurbsFactory::create_circle(1., nc);

	Point3 direction(0, 0, 3);
	NurbsSurface ns;

	// extrude
	NurbsExtrude ne;
	ne.extrude(nc, direction, ns);

	Mesh m;
	NurbsUtil::to_mesh(ns,m,10);
	OBJFile::save("test_nurbsextrude_cylinder.obj", m);

	StepWriter sw;
	sw.open("test_nurbsextrude_cylinder.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsextrude_random_deg2()
{
	cout << endl << "test_nurbsextrude_random_deg2" << endl;

	NurbsCurve nc;
	int degree = 2;
	vector<Point3> points;
	int nbPoints = 5;

	for (int i = 0; i < nbPoints; i++)
		points.push_back(Point3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX));

	NurbsUtil::create_curve_from_points(points, degree,nc);

	Point3 direction(0, 0, 1.);
	NurbsSurface ns;

	// extrude
	NurbsExtrude ne;
	ne.extrude(nc, direction, ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 10);
	OBJFile::save("test_nurbsextrude_random_deg2.obj", m);

	StepWriter sw;
	sw.open("test_nurbsextrude_random_deg2.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsextrude_face()
{
	cout << endl << "test_nurbsextrude_face" << endl;

	// Create a NURBS cube
	NurbsSolid cube;
	NurbsFactory::create_box(10., 10., 10., cube);

	int originalSurfaceCount = cube.surfaces().size();

	// Compute extrusion direction using face normal at center
	const auto& surfaces = cube.surfaces();
	Point3 normal;
	if (!surfaces[0].normal(0.5, 0.5, normal)) {
		cout << "Failed to compute normal" << endl;
		return;
	}
	Point3 direction = normal * 5.0;

	// Extrude the top face (index 0)
	NurbsExtrude ne;
	bool success = ne.extrude_face(cube, 0, direction);

	if (!success) {
		cout << "Failed to extrude face" << endl;
		return;
	}

	if (cube.surfaces().size() != originalSurfaceCount + 1) {
		cout << "Surface count not increased correctly" << endl;
		return;
	}

	cout << "Face extrusion test passed" << endl;

	// Save the result
	Mesh m;
	NurbsUtil::to_mesh(cube, m, 16);
	OBJFile::save("test_nurbsextrude_face.obj", m);

	StepWriter sw;
	sw.open("test_nurbsextrude_face.step");
	sw.write(cube);
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	test_nurbsextrude_cylinder();
	test_nurbsextrude_random_deg2();
	test_nurbsextrude_face();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////