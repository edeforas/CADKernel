#include "NurbsExtrude.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
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
int main()
{
	test_nurbsextrude_cylinder();
	test_nurbsextrude_random_deg2();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////