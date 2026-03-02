#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsExtrude.h"
#include "NurbsUtil.h"
#include "Transform.h"

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
	sw.open(filename);
	test_bool(sw.is_open(), string("step writer should open file: ") + filename);
	sw.write(n);
	sw.close();

	ifstream in(filename.c_str());
	test_bool(in.is_open(), string("step file should exist: ") + filename);

	OBJWriter o;
	Mesh m;
	NurbsUtil::to_mesh(n, m);
	o.open(filename + ".obj");
	o.write(m);
}

///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_circle()
{
	cout << endl << "test_nurbsfactory_create_circle" << endl;

	NurbsCurve n;

	for (double dRadius = 0.1; dRadius < 100.; dRadius *= 2.)
	{
		NurbsFactory::create_circle(dRadius, n);

		for (double u = 0.; u <= 1.; u += 0.01)
		{
			Point3 p;
			n.evaluate(u, p);
			double norm = p.norm();
			//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
			test_near(norm, dRadius, 1.e-10);
		}
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_sphere()
{
	cout << endl << "test_nurbsfactory_create_sphere" << endl;

	for (double dRadius = 0.1; dRadius < 100.; dRadius *= 2.)
	{
		NurbsSurface ns;
		NurbsFactory::create_sphere(dRadius, ns);

		for (double u = 0.; u <= 1.; u += 0.01)
			for (double v = 0.; v <= 1.; v += 0.01)
			{
				Point3 p;
				ns.evaluate(u, v, p);
				double norm = p.norm();
				//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
				test_near(norm, dRadius, 1.e-10);
			}
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_solid_cylinder()
{
	cout << endl << "test_nurbsfactory_create_solid_cylinder" << endl;

	NurbsSolid ns;
	NurbsFactory::create_cylinder(10, 30, ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 10);
	OBJFile::save("test_nurbsfactory_create_solid_cylinder.obj", m);

	save_solid(ns, "test_nurbsfactory_create_solid_cylinder");
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_solid_sphere()
{
	cout << endl << "test_nurbsfactory_create_solid_sphere" << endl;

	NurbsSolid n;
	NurbsFactory::create_sphere(20, n);

	save_solid(n, "test_nurbsfactory_create_solid_sphere");
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_solid_torus()
{
	cout << endl << "test_nurbsfactory_create_solid_torus" << endl;

	NurbsSolid n;
	NurbsFactory::create_torus(30,10, n);

	save_solid(n, "test_nurbsfactory_create_solid_torus");
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_solid_box()
{
	cout << endl << "test_nurbsfactory_create_solid_box" << endl;

	NurbsSolid ns;
	NurbsFactory::create_box(10,40,90, ns);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 3);
	OBJFile::save("test_nurbsfactory_create_solid_box.obj", m);

	StepWriter sw;
	sw.open("test_nurbsfactory_create_solid_box.step");
	sw.write(ns);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_fill_curve_to_surface_export()
{
	cout << endl << "test_nurbsfactory_fill_curve_to_surface_export" << endl;

	NurbsCurve nc;
	NurbsFactory::create_circle(25., nc);

	NurbsSurface ns;
	NurbsFactory::create_filled_surface(nc, ns);

	test_bool(ns.nb_points_u() == nc.nb_points(), "filled surface should preserve number of U control points");
	test_bool(ns.nb_points_v() == 2, "filled surface should have two V rows");

	const string stepFilename = "test_nurbsfactory_fill_curve_to_surface.step";
	StepWriter sw;
	sw.open(stepFilename);
	test_bool(sw.is_open(), string("step writer should open file: ") + stepFilename);
	sw.write(ns);
	sw.close();

	ifstream in(stepFilename.c_str());
	test_bool(in.is_open(), string("step file should exist: ") + stepFilename);

	Mesh m;
	NurbsUtil::to_mesh(ns, m, 30);
	test_bool(OBJFile::save("test_nurbsfactory_fill_curve_to_surface.obj", m), "obj export should succeed");
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsfactory_create_solid_torus_transform()
{
	cout << endl << "test_nurbsfactory_create_solid_torus_transform" << endl;

	OBJWriter ow;
	ow.open("test_nurbsfactory_create_solid_torus_transform.obj");

	StepWriter sw;
	sw.open("test_nurbsfactory_create_solid_torus_transform.step");

	for (double i=0.5; i < 5; i *= 2.)
		for (double j = 0.5; j < 5; j *= 2.)
			for (double k = 0.5; k < 5; k *= 2.)
			{
				NurbsSurface ns;
				NurbsFactory::create_torus(30, 10, ns);
				Scale scale(i, j, k);
				Translation translate(Point3(i, j, k) * 100.);
				RotationEulerAngles rotate(i*10, j*10, k*10);

				rotate.apply_all(ns.points());
				scale.apply_all(ns.points());
				translate.apply_all(ns.points());

				Mesh m;
				NurbsUtil::to_mesh(ns, m, 4);
				ow.write(m);
				sw.write(ns);
		}

	ow.close();
	sw.close();
}
///////////////////////////////////////////////////////////////////////////

int main()
{
	test_nurbsfactory_create_circle();
	test_nurbsfactory_create_sphere();
	test_nurbsfactory_create_solid_cylinder();
	test_nurbsfactory_create_solid_sphere();
	test_nurbsfactory_create_solid_torus();
	test_nurbsfactory_create_solid_box();
	test_nurbsfactory_fill_curve_to_surface_export();
	test_nurbsfactory_create_solid_torus_transform();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////