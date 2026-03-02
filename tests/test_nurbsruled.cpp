#include "NurbsRuled.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "NurbsFactory.h"
#include "Transform.h"

#include "OBJFile.h"
#include "StepFile.h"

#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

void test_near(double a, double ref, double epsilon=1.e-10,const string& sMessage="")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsruled_deg1()
{
	cout << endl << "test_nurbsruled_deg1" << endl;

	// test handcrafted ruled surface
	NurbsSurface n;
	int degree = 1;
	vector<Point3> points = { Point3(0.,0.,1.),Point3(1.,0.,-1.),Point3(0.,1.,-1.),Point3(1.,1.,1.) };

	n.set_degree(degree, degree);
	n.set_points(points, 2, 2);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 p;
			n.evaluate(u, v, p);
			test_near(p.x(), u, 1e-10);
			test_near(p.y(), v, 1e-10);

			test_bool(p.z() <= 1.);
			test_bool(p.z() >= -1.);

			//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << endl;
		}

	Mesh m;
	NurbsUtil::to_mesh(n, m, 10);
	OBJFile::save("test_nurbsruled_deg1.obj", m);

	StepWriter sw;
	sw.open("test_nurbsruled_deg1.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsruled_deg2()
{
	cout << endl << "test_nurbsruled_deg2" << endl;

	// test handcrafted ruled surface
	NurbsSurface n;
	int degreeU = 2;
	vector<Point3> points = {
		Point3(0.,0.,0.),Point3(1.,1.,0.),Point3(2.,0.,0.),
		Point3(0.,1.,-1.),Point3(1.,0.,-1.),Point3(2.,1.,-1.)
	};

	n.set_degree(1, 1);
	n.set_points(points, 3, 2);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();

	Mesh m;
	NurbsUtil::to_mesh(n, m, 10);
	OBJFile::save("test_nurbsruled_deg2.obj", m);

	StepWriter sw;
	sw.open("test_nurbsruled_deg_etienne.step");
	sw.write(n);

}
///////////////////////////////////////////////////////////////////////////
void test_nurbsruled_cylinder()
{
	cout << endl << "test_nurbsruled_cylinder" << endl;

	//create profile curve
	NurbsCurve nc1,nc2;
	NurbsFactory::create_circle(1., nc1);
	NurbsFactory::create_circle(1., nc2);

	// translate nc2
	Translation t(Point3(0, 0, 1));
	t.apply_all(nc2.points());

	NurbsSurface n;

	// ruled
	NurbsRuled ne;
	ne.create_ruled(nc1,nc2, n);

	Mesh m;
	NurbsUtil::to_mesh(n,m,10);
	OBJFile::save("test_nurbsruled_cylinder.obj", m);

	StepWriter sw;
	sw.open("test_nurbsruled_cylinder.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsruled_ribbon_deg3()
{
	NurbsCurve n1, n2;
	vector<Point3> points1, points2;
	int nbPoints = 11;
	for (int i = 0; i < nbPoints; i++)
		points1.push_back(Point3((double)rand() / RAND_MAX + i * 2, (double)rand() / RAND_MAX, 10.*(double)rand() / RAND_MAX));
	for (int i = 0; i < nbPoints; i++)
		points2.push_back(Point3((double)rand() / RAND_MAX + i * 2, (double)rand() / RAND_MAX + 2, 10.*(double)rand() / RAND_MAX));

	int degree = 3;
	NurbsUtil::create_curve_from_points(points1, degree, n1);
	NurbsUtil::create_curve_from_points(points2, degree, n2);

	NurbsSurface n;

	// ruled
	NurbsRuled ne;
	ne.create_ruled(n1, n2, n);

	Mesh m;
	NurbsUtil::to_mesh(n, m, 10);
	OBJFile::save("test_nurbsruled_ribbon_deg3.obj", m);

	StepWriter sw;
	sw.open("test_nurbsruled_ribbon_deg3.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
int main()
{

	test_nurbsruled_deg2();


	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////