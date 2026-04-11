#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSweep.h"
#include "NurbsFactory.h"
#include "NurbsUtil.h"

#include <iostream>
#include <cmath>
#include <limits>
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

void test_near(double a, double ref, double epsilon = 1.e-10, const string& sMessage = "")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << " value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}

void test_sweep_linear_segment_along_line()
{
	cout << endl << "test_sweep_linear_segment_along_line" << endl;

	NurbsCurve profile;
	profile.set_degree(1);
	profile.set_points({ Point3(0., 0., 0.), Point3(2., 0., 0.) });
	profile.set_knots({ 0., 0., 1., 1. });
	profile.set_weights({ 1., 1. });

	NurbsCurve path;
	path.set_degree(1);
	path.set_points({ Point3(0., 0., 0.), Point3(0., 0., 3.) });
	path.set_knots({ 0., 0., 1., 1. });
	path.set_weights({ 1., 1. });

	NurbsSurface surface;
	bool ok = NurbsSweep::sweep(profile, path, surface);
	test_bool(ok, "sweep should succeed for valid linear curves");

	for (double u = 0.; u <= 1.; u += 0.2)
		for (double v = 0.; v <= 1.; v += 0.2)
		{
			Point3 p;
			surface.evaluate(u, v, p);
			test_near(p.x(), 2. * u, 1.e-10, "x mismatch");
			test_near(p.y(), 0., 1.e-10, "y mismatch");
			test_near(p.z(), 3. * v, 1.e-10, "z mismatch");
		}
}

void test_sweep_recovers_from_invalid_knots_weights_and_degree()
{
	cout << endl << "test_sweep_recovers_from_invalid_knots_weights_and_degree" << endl;

	NurbsCurve profile;
	profile.set_degree(5); // intentionally invalid for 2 points
	profile.set_points({ Point3(0., 1., 0.), Point3(1., 1., 0.) });
	profile.set_knots({ 1., 0. }); // intentionally invalid
	profile.set_weights({ std::numeric_limits<double>::quiet_NaN(), 0.0 }); // intentionally invalid

	NurbsCurve path;
	path.set_degree(4); // intentionally invalid for 2 points
	path.set_points({ Point3(0., 0., 0.), Point3(0., 2., 0.) });
	path.set_knots({ 1., 0. }); // intentionally invalid
	path.set_weights({ -2.0, std::numeric_limits<double>::quiet_NaN() }); // intentionally invalid

	NurbsSurface surface;
	bool ok = NurbsSweep::sweep(profile, path, surface);
	test_bool(ok, "sweep should recover from malformed curve data");

	Point3 p00, p11;
	surface.evaluate(0., 0., p00);
	surface.evaluate(1., 1., p11);

	test_bool(std::isfinite(p00.x()) && std::isfinite(p00.y()) && std::isfinite(p00.z()), "p00 should be finite");
	test_bool(std::isfinite(p11.x()) && std::isfinite(p11.y()) && std::isfinite(p11.z()), "p11 should be finite");

	test_near(p00.x(), 0., 1.e-10, "p00.x mismatch");
	test_near(p00.y(), 1., 1.e-10, "p00.y mismatch");
	test_near(p00.z(), 0., 1.e-10, "p00.z mismatch");

	test_near(p11.x(), 1., 1.e-10, "p11.x mismatch");
	test_near(p11.y(), 3., 1.e-10, "p11.y mismatch");
	test_near(p11.z(), 0., 1.e-10, "p11.z mismatch");
}

void test_sweep_propagates_closed_flags()
{
	cout << endl << "test_sweep_propagates_closed_flags" << endl;

	NurbsCurve profile;
	profile.set_degree(1);
	profile.set_points({ Point3(0., 0., 0.), Point3(1., 0., 0.), Point3(0., 0., 0.) });
	profile.set_uniform();
	profile.set_equals_weights();

	NurbsCurve path;
	path.set_degree(1);
	path.set_points({ Point3(0., 0., 0.), Point3(0., 1., 0.), Point3(0., 0., 0.) });
	path.set_uniform();
	path.set_equals_weights();

	NurbsSurface surface;
	bool ok = NurbsSweep::sweep(profile, path, surface);
	test_bool(ok, "sweep should succeed for closed profile/path");

	test_bool(surface.is_closed_u(), "surface closed_u should follow profile closure");
	test_bool(surface.is_closed_v(), "surface closed_v should follow path closure");

	Point3 p;
	surface.evaluate(0.3, 0.7, p);
	test_bool(std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()), "evaluated point should be finite");
}

void test_sweep_perpendicular_rotates_profile_plane()
{
	cout << endl << "test_sweep_perpendicular_rotates_profile_plane" << endl;

	NurbsCurve profile;
	profile.set_degree(1);
	profile.set_points({ Point3(1., 1., 0.), Point3(-1., 1., 0.) });
	profile.set_knots({ 0., 0., 1., 1. });
	profile.set_weights({ 1., 1. });

	NurbsCurve path;
	path.set_degree(1);
	path.set_points({ Point3(0., 0., 0.), Point3(0., 1., 0.) });
	path.set_knots({ 0., 0., 1., 1. });
	path.set_weights({ 1., 1. });

	NurbsSurface surface;
	bool ok = NurbsSweep::sweep(profile, path, surface, true);
	test_bool(ok, "perpendicular sweep should succeed");

	Point3 p;
	surface.evaluate(0.25, 0.5, p);
	test_bool(std::abs(p.z()) > 1.e-6, "perpendicular sweep should produce non-planar points");
}

void test_sweep_perpendicular_minimizes_torsion()
{
	cout << endl << "test_sweep_perpendicular_minimizes_torsion" << endl;

	NurbsCurve profile;
	profile.set_degree(1);
	profile.set_points({ Point3(1., 0., 0.), Point3(-1., 0., 0.) });
	profile.set_knots({ 0., 0., 1., 1. });
	profile.set_weights({ 1., 1. });

	NurbsCurve path;
	path.set_degree(1);
	path.set_points({
		Point3(1., 0., 0.),
		Point3(0.7, 0.7, 0.5),
		Point3(0., 1., 1.),
		Point3(-0.7, 0.7, 1.5),
		Point3(-1., 0., 2.)
	});
	path.set_uniform();
	path.set_equals_weights();

	NurbsSurface surface;
	bool ok = NurbsSweep::sweep(profile, path, surface, true);
	test_bool(ok, "perpendicular sweep should succeed");

	Point3 prevDir;
	bool hasPrev = false;
	for (int j = 0; j < 5; ++j) {
		double v = j / 4.0;
		Point3 p0, p1;
		surface.evaluate(0.0, v, p0);
		surface.evaluate(1.0, v, p1);
		Point3 dir = p1 - p0;
		dir.normalize();
		if (hasPrev) {
			double alignment = dir.dot_product(prevDir);
			test_bool(alignment > 0.8, "consecutive profile orientations should remain smooth");
		}
		prevDir = dir;
		hasPrev = true;
	}
}

int main()
{
	test_sweep_linear_segment_along_line();
	test_sweep_recovers_from_invalid_knots_weights_and_degree();
	test_sweep_propagates_closed_flags();
	test_sweep_perpendicular_rotates_profile_plane();
	test_sweep_perpendicular_minimizes_torsion();

	cout << "Test Finished.";
	return 0;
}
