#include "NurbsRevolve.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsFactory.h"
#include "NurbsUtil.h"

#include "MeshUtil.h"
#include "OBJFile.h"
#include "StepFile.h"

#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

#ifndef M_PI
#define M_PI (3.1415926535897932384)
#endif
void test_near(double a, double ref, double epsilon=1.e-10,const string& sMessage="")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_flatdisk()
{
	cout << endl << "test_nurbsrevolve_flatdisk" << endl;

	//create profile curve
	vector<Point3> points = {Point3(0.,0.,0.),Point3(1.,0.,0.) };
	NurbsCurve nc;
	NurbsUtil::create_curve_from_points(points, 1, nc);

	// revolve
	NurbsRevolve nr;
	NurbsSurface n;
	nr.revolve(nc, n);

	Mesh m;
	NurbsUtil::to_mesh(n, m, 10);

	OBJWriter ow;
	ow.open("test_nurbsrevolve_flatdisk.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_flatdisk.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_cylinder()
{
	cout << endl << "test_nurbsrevolve_cylinder" << endl;

	//create profile curve
	vector<Point3> points = { Point3(0.,0.,0.),Point3(1.,0.,0.),Point3(1.,0.,3.),Point3(0.,0.,3.) };
	NurbsCurve nc;
	NurbsUtil::create_curve_from_points(points, 1, nc);

	// revolve
	NurbsRevolve nr;
	NurbsSurface n;
	nr.revolve(nc, n);

	Mesh m;
	NurbsUtil::to_mesh(n, m, 10);

	OBJWriter ow;
	ow.open("test_nurbsrevolve_cylinder.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_cylinder.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_vase()
{
	cout << endl << "test_nurbsrevolve_vase" << endl;

	//create profile curve
	vector<Point3> points = {
		Point3(0.,0.,0.),Point3(1.,0.,0.),
		Point3(1.,0.,0.),Point3(2.,1.,1.),Point3(1.,1.,2.),
		Point3(2.,2.,3.),Point3(3.,0.,4.),Point3(0.,2.,5.)
	};
	NurbsCurve nc;
	NurbsUtil::create_curve_from_points(points, 3, nc);

	// revolve
	NurbsRevolve nr;
	NurbsSurface n;
	nr.revolve(nc, n);

	Mesh m;
	NurbsUtil::to_mesh(n,m,10);

	OBJWriter ow;
	ow.open("test_nurbsrevolve_vase.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_vase.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_sphere()
{
	cout << endl << "test_nurbsrevolve_sphere" << endl;

	//creation with elementary operators
	//create profile curve: half circle
	NurbsCurve nc;
	vector<Point3> points = {
		Point3(0., 0.,-1.),Point3(0.,1.,-1.),Point3(0., 1.,0.),Point3(0.,1.,1.),Point3(0.,0.,1.)
	};
	vector<double> weights = { 1. , 1. / sqrt(2.),1., 1. / sqrt(2.), 1. };
	vector<double> knots = { 0., 0., 0., 0.5, 0.5, 1., 1., 1. };
	nc.set_degree(2);
	nc.set_points(points);
	nc.set_weights(weights);
	nc.set_knots(knots);

	NurbsRevolve nr;
	NurbsSurface ns;
	nr.revolve(nc, ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 10);

	// test all points are on distance 1. to center
	for (int i = 0; i < m.nb_vertices(); i++)
	{
		Point3 p;
		m.get_vertex(i, p);
		test_near(p.norm(), 1.);
	}

	OBJWriter ow;
	ow.open("test_nurbsrevolve_sphere.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_sphere.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_create_sphere()
{
	// use nurbs factory that call revolve
	cout << endl << "test_nurbsrevolve_create_sphere" << endl;

	NurbsSurface ns;
	NurbsFactory::create_sphere(1., ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 10);

	// test all points are on distance 1. to center
	for (int i = 0; i < m.nb_vertices(); i++)
	{
		Point3 p;
		m.get_vertex(i, p);
		test_near(p.norm(), 1.);
	}

	OBJWriter ow;
	ow.open("test_nurbsrevolve_create_sphere.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_create_sphere.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_create_torus()
{
	// call nurbs factory
	cout << endl << "test_nurbsrevolve_create_torus" << endl;

	NurbsSurface ns;
	NurbsFactory::create_torus(5., 2., ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 10);

	OBJWriter ow;
	ow.open("test_nurbsrevolve_create_torus.obj");
	ow.write(m);

	StepWriter sw;
	sw.open("test_nurbsrevolve_create_torus.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsrevolve_any_angle_step_export()
{
	cout << endl << "test_nurbsrevolve_any_angle_step_export" << endl;

	vector<Point3> points = { Point3(1.,0.,-1.), Point3(1.,0.,1.) };
	NurbsCurve nc;
	NurbsUtil::create_curve_from_points(points, 1, nc);

	NurbsRevolve nr;
	const vector<double> angles = {
		M_PI / 6.,
		M_PI,
		1.75 * M_PI,
		2. * M_PI
	};

	for (int i = 0; i < (int)angles.size(); ++i)
	{
		NurbsSurface ns;
		bool ok = nr.revolve(nc, angles[i], ns);
		assert(ok);

		test_near(std::fabs(ns.knots_u().front()), 0., 1.e-12, "u knots should start at 0");
		test_near(std::fabs(ns.knots_u().back() - 1.), 0., 1.e-12, "u knots should end at 1");

		if (angles[i] >= 2. * M_PI - 1.e-9)
			assert(ns.is_closed_u());
		else
			assert(!ns.is_closed_u());

		StepWriter sw;
		string filename = string("test_nurbsrevolve_angle_") + to_string(i) + string(".step");
		sw.open(filename);
		assert(sw.is_open());
		sw.write(ns);
		sw.close();

		OBJWriter o;
		Mesh m;
		NurbsUtil::to_mesh(ns, m);
		o.open(filename + ".obj");
		o.write(m);

	}

	NurbsSurface nsDeg;
	bool okDeg = nr.revolve_deg(nc, 270., nsDeg);
	assert(okDeg);
	assert(!nsDeg.is_closed_u());

	StepWriter swDeg;
	swDeg.open("test_nurbsrevolve_angle_deg_270.step");
	assert(swDeg.is_open());
	swDeg.write(nsDeg);
	swDeg.close();
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	test_nurbsrevolve_cylinder();
	test_nurbsrevolve_flatdisk();

	test_nurbsrevolve_vase();
	test_nurbsrevolve_sphere();
	test_nurbsrevolve_create_sphere();
	test_nurbsrevolve_create_torus();
	test_nurbsrevolve_any_angle_step_export();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////