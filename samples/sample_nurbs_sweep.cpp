#include "NurbsSweep.h"
#include "NurbsFactory.h"
#include "NurbsCurve.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "OBJFile.h"
#include "Mesh.h"

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main()
{
	cout << "NURBS Sweep Sample with Endpoint Closures" << endl;
	cout << "This sample demonstrates sweeping a circular profile along a helical path." << endl;
	cout << "It creates sweeps with different endpoint closure options." << endl;
	cout << endl;

	// Create a profile curve: a circle
	NurbsCurve profile;
	NurbsFactory::create_circle(1.0, profile); // radius 1
	cout << "Created circular profile with radius 1.0" << endl;

	// Create a path curve: a helical 3D curve
	NurbsCurve path;
	std::vector<Point3> pathPoints;
	int numPoints = 20;
	for (int i = 0; i < numPoints; ++i) {
		double t = (double)i / (numPoints - 1);
		double angle = t * 2.0 * acos(-1.0) * 2.0; // two full turns
		double radius = 3.0;
		double x = radius * cos(angle);
		double y = radius * sin(angle);
		double z = t * 10.0;
		pathPoints.push_back(Point3(x, y, z));
	}
	NurbsUtil::create_curve_from_points(pathPoints, 3, path); // degree 3 for smooth curve
	cout << "Created helical path with " << numPoints << " control points" << endl;
	cout << endl;

	// Sweep with mixed caps (both orientations)
	cout << "Creating sweeps with mixed caps..." << endl;

	NurbsSolid solidPerp;
	NurbsSweep::sweep_solid(profile, path, solidPerp, true,
		NurbsSweep::EndpointClosure::HalfSphere,
		NurbsSweep::EndpointClosure::Disk);
	
	OBJWriter ow;
	Mesh meshPerp;
	NurbsUtil::to_mesh(solidPerp, meshPerp, 5);
	ow.open("sample_nurbs_sweep_mixed_caps.obj");
	ow.write(meshPerp);
	ow.close();
	
	StepWriter sw;
	sw.open("sample_nurbs_sweep_mixed_caps.step");
	sw.write(solidPerp);
	sw.close();

	return 0;
}