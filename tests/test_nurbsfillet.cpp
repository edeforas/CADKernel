#include "NurbsFillet.h"

#include "NurbsCurve.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "StepFile.h"

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

bool file_exists(const string& filename)
{
	ifstream in(filename.c_str());
	return in.is_open();
}

void test_nurbs_chamfer_surface_and_step()
{
	cout << endl << "test_nurbs_chamfer_surface_and_step" << endl;

	NurbsCurve c1, c2;
	NurbsUtil::create_curve_from_points({ Point3(0., 0., 0.), Point3(2., 0., 1.), Point3(4., 0., 0.) }, 2, c1);
	NurbsUtil::create_curve_from_points({ Point3(0., 2., 0.), Point3(2., 2., 1.), Point3(4., 2., 0.) }, 2, c2);

	NurbsFillet nf;
	NurbsSurface chamfer;
	bool ok = nf.create_chamfer(c1, c2, 0.2, 0.3, chamfer);
	if (!ok)
		ok = nf.create_chamfer(c1, c2, 0.2, 0.3, chamfer);
	test_bool(ok, "chamfer should be created");

	test_bool(chamfer.degree_v() == 1, "chamfer should be linear in v");
	test_bool(chamfer.nb_points_v() == 2, "chamfer should have 2 rows in v");
	test_bool(chamfer.nb_points_u() == c1.nb_points(), "chamfer should preserve u control points");

	StepWriter sw;
	sw.open("test_nurbs_chamfer.step");
	test_bool(sw.is_open(), "step writer should open for chamfer");
	sw.write(chamfer);
	sw.close();
	test_bool(file_exists("test_nurbs_chamfer.step"), "chamfer step file should exist");
}

void test_nurbs_fillet_surface_and_step()
{
	cout << endl << "test_nurbs_fillet_surface_and_step" << endl;

	NurbsCurve c1, c2;
	NurbsUtil::create_curve_from_points({ Point3(0., 0., 0.), Point3(2., 0., 1.), Point3(4., 0., 0.) }, 2, c1);
	NurbsUtil::create_curve_from_points({ Point3(0., 2., 0.), Point3(2., 2., 1.), Point3(4., 2., 0.) }, 2, c2);

	NurbsFillet nf;
	NurbsSurface fillet;
	bool ok = nf.create_fillet(c1, c2, 0.25, fillet);
	test_bool(ok, "fillet should be created");

	test_bool(fillet.degree_v() == 2, "fillet should be quadratic in v");
	test_bool(fillet.nb_points_v() == 3, "fillet should have 3 rows in v");
	test_bool(fillet.nb_points_u() == c1.nb_points(), "fillet should preserve u control points");

	StepWriter sw;
	sw.open("test_nurbs_fillet.step");
	test_bool(sw.is_open(), "step writer should open for fillet");
	sw.write(fillet);
	sw.close();
	test_bool(file_exists("test_nurbs_fillet.step"), "fillet step file should exist");
}

void test_nurbs_fillet_incompatible_curves()
{
	cout << endl << "test_nurbs_fillet_incompatible_curves" << endl;

	NurbsCurve c1, c2;
	NurbsUtil::create_curve_from_points({ Point3(0., 0., 0.), Point3(1., 0., 0.), Point3(2., 0., 0.) }, 2, c1);
	NurbsUtil::create_curve_from_points({ Point3(0., 1., 0.), Point3(2., 1., 0.) }, 1, c2);

	NurbsFillet nf;
	NurbsSurface out;
	test_bool(!nf.create_chamfer(c1, c2, 0.1, 0.1, out), "incompatible chamfer input should fail");
	test_bool(!nf.create_fillet(c1, c2, 0.1, out), "incompatible fillet input should fail");
}

void test_nurbs_fillet_chamfer_on_solid_edges_step()
{
	cout << endl << "test_nurbs_fillet_chamfer_on_solid_edges_step" << endl;

	NurbsSolid solid;
	NurbsFactory::create_cylinder(2., 4., solid);

	NurbsFillet nf;
	NurbsSurface fillet, chamfer;

	bool okFillet = nf.create_fillet_on_solid(
		solid,
		0, NurbsFillet::EdgeVMax,
		1, NurbsFillet::EdgeVMin,
		0.2,
		fillet);
	test_bool(okFillet, "solid edge fillet should be created");
	test_bool(fillet.nb_points_u() > 0 && fillet.nb_points_v() == 3, "solid fillet topology should be valid");

	bool okChamfer = nf.create_chamfer_on_solid(
		solid,
		0, NurbsFillet::EdgeVMax,
		1, NurbsFillet::EdgeVMin,
		0.1, 0.1,
		chamfer);
	test_bool(okChamfer, "solid edge chamfer should be created");
	test_bool(chamfer.nb_points_u() > 0 && chamfer.nb_points_v() == 2, "solid chamfer topology should be valid");

	StepWriter swF;
	swF.open("test_nurbs_solid_edge_fillet.step");
	test_bool(swF.is_open(), "step writer should open for solid edge fillet");
	swF.write(fillet);
	swF.close();
	test_bool(file_exists("test_nurbs_solid_edge_fillet.step"), "solid edge fillet step file should exist");

	StepWriter swC;
	swC.open("test_nurbs_solid_edge_chamfer.step");
	test_bool(swC.is_open(), "step writer should open for solid edge chamfer");
	swC.write(chamfer);
	swC.close();
	test_bool(file_exists("test_nurbs_solid_edge_chamfer.step"), "solid edge chamfer step file should exist");
}

void test_nurbs_fillet_chamfer_on_solid_edges_auto_step()
{
	cout << endl << "test_nurbs_fillet_chamfer_on_solid_edges_auto_step" << endl;

	NurbsSolid solid;
	NurbsFactory::create_cylinder(2., 4., solid);

	NurbsFillet nf;
	NurbsSurface fillet, chamfer;

	bool okFillet = nf.create_fillet_on_solid_auto(solid, 0, 1, 0.2, fillet);
	test_bool(okFillet, "auto solid edge fillet should be created");
	test_bool(fillet.nb_points_u() > 0 && fillet.nb_points_v() == 3, "auto solid fillet topology should be valid");

	bool okChamfer = nf.create_chamfer_on_solid_auto(solid, 0, 1, 0.1, 0.1, chamfer);
	test_bool(okChamfer, "auto solid edge chamfer should be created");
	test_bool(chamfer.nb_points_u() > 0 && chamfer.nb_points_v() == 2, "auto solid chamfer topology should be valid");

	StepWriter swF;
	swF.open("test_nurbs_solid_edge_fillet_auto.step");
	test_bool(swF.is_open(), "step writer should open for auto solid edge fillet");
	swF.write(fillet);
	swF.close();
	test_bool(file_exists("test_nurbs_solid_edge_fillet_auto.step"), "auto solid edge fillet step file should exist");

	StepWriter swC;
	swC.open("test_nurbs_solid_edge_chamfer_auto.step");
	test_bool(swC.is_open(), "step writer should open for auto solid edge chamfer");
	swC.write(chamfer);
	swC.close();
	test_bool(file_exists("test_nurbs_solid_edge_chamfer_auto.step"), "auto solid edge chamfer step file should exist");
}

int main()
{
	test_nurbs_chamfer_surface_and_step();
	test_nurbs_fillet_surface_and_step();
	test_nurbs_fillet_incompatible_curves();
	test_nurbs_fillet_chamfer_on_solid_edges_step();
	test_nurbs_fillet_chamfer_on_solid_edges_auto_step();

	cout << "Test Finished.";
	return 0;
}
