#include "NurbsCurve.h"
#include "NurbsCurveJoin.h"
#include "NurbsFactory.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "OBJFile.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

void sample_curve_to_polyline(const NurbsCurve& c, std::vector<Point3>& polyline, int samples = 60)
{
	polyline.clear();
	if (samples < 2)
		samples = 2;
	for (int i = 0; i <= samples; ++i)
	{
		double u = (double)i / samples;
		Point3 p;
		c.evaluate(u, p);
		polyline.push_back(p);
	}
}

bool createConnectorCurve(const NurbsCurve& c1, const NurbsCurve& c2, int continuity, NurbsCurve& out)
{
	if (c1.nb_points() < 2 || c2.nb_points() < 2)
		return false;

	Point3 p0 = c1.points().back();
	Point3 p3 = c2.points().front();
	out.clear();

	if (continuity == 0)
	{
		std::vector<Point3> pts = { p0, p3 };
		out.set_degree(1);
		out.set_points(pts);
		out.set_equals_weights();
		out.set_knots({ 0., 0., 1., 1. });
		return true;
	}

	Point3 t0, t1;
	if (!c1.tangent(1.0, t0))
		return false;
	if (!c2.tangent(0.0, t1))
		return false;

	if (t0.norm_square() < 1e-12 || t1.norm_square() < 1e-12)
		return false;

	double d = std::sqrt(p0.distance_square(p3));
	double h = std::max(1e-6, d * 0.4);

	if (continuity == 1)
	{
		Point3 q1 = p0 + t0 * h;
		Point3 q2 = p3 - t1 * h;
		std::vector<Point3> pts = { p0, q1, q2, p3 };
		out.set_degree(3);
		out.set_points(pts);
		out.set_equals_weights();
		out.set_knots({ 0.,0.,0.,0.,1.,1.,1.,1. });
		return true;
	}

	if (continuity == 2)
	{
		// G2: use a degree 4 connecting curve with extra control point for smooth curvature.
		Point3 q1 = p0 + t0 * (h * 0.5);
		Point3 q3 = p3 - t1 * (h * 0.5);
		Point3 q2 = (p0 + p3) * 0.5 + (t0 - t1) * (d * 0.05);
		std::vector<Point3> pts = { p0, q1, q2, q3, p3 };
		out.set_degree(4);
		out.set_points(pts);
		out.set_equals_weights();
		out.set_knots({ 0.,0.,0.,0.,0.,1.,1.,1.,1.,1. });
		return true;
	}

	return false;
}

int main(int argc, char* argv[])
{
	int continuity = 1; // default G1
	if (argc > 1)
	{
		string arg = argv[1];
		if (arg == "G0" || arg == "g0") continuity = 0;
		else if (arg == "G1" || arg == "g1") continuity = 1;
		else if (arg == "G2" || arg == "g2") continuity = 2;
		else
		{
			cout << "Usage: " << argv[0] << " [G0|G1|G2]" << endl;
			return 1;
		}
	}

	cout << "sample_nurbs_curve_continuity: connecting two curves with continuity G" << continuity << endl;

	// create two sample curves (non colinear) to connect
	NurbsCurve c1, c2;
	NurbsUtil::create_curve_from_points({ Point3(0, 0, 0), Point3(2, 3, 0), Point3(4, 3, 0) }, 3, c1);
	NurbsUtil::create_curve_from_points({ Point3(8, 0, 0), Point3(10, -3, 0), Point3(12, -1, 0) }, 3, c2);

	NurbsCurve connector;
	if (!NurbsCurveJoin::create_connector(c1, c2, (NurbsCurveJoin::Continuity)continuity, connector))
	{
		cout << "Failed to create connecting curve" << endl;
		return 1;
	}

	NurbsCurve chamferCurve, filletCurve;
	if (!NurbsCurveJoin::create_chamfer(c1, c2, 1.0, chamferCurve))
	{
		cout << "Failed to create chamfer curve" << endl;
		return 1;
	}
	if (!NurbsCurveJoin::create_fillet(c1, c2, 1.0, filletCurve))
	{
		cout << "Failed to create fillet curve" << endl;
		return 1;
	}

	// Save curves to STEP
	StepWriter sw;
	sw.open("sample_nurbs_curve_continuity.step");
	sw.write(c1);
	sw.write(c2);
	sw.write(connector);
	sw.write(chamferCurve);
	sw.write(filletCurve);
	sw.close();

	// Write to OBJ as polylines
	OBJWriter ow;
	ow.open("sample_nurbs_curve_continuity.obj");
	vector<Point3> pl;
	sample_curve_to_polyline(c1, pl, 60);
	ow.write(pl);
	sample_curve_to_polyline(c2, pl, 60);
	ow.write(pl);
	sample_curve_to_polyline(connector, pl, 80);
	ow.write(pl);
	sample_curve_to_polyline(chamferCurve, pl, 2);
	ow.write(pl);
	sample_curve_to_polyline(filletCurve, pl, 80);
	ow.write(pl);
	ow.close();

	cout << "Wrote sample_nurbs_curve_continuity.step and sample_nurbs_curve_continuity.obj" << endl;
	cout << "(use command line argument G0, G1 or G2)" << endl;
	return 0;
}
