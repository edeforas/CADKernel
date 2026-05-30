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
	cout << "NURBS Sweep Sample" << endl;
	cout << "Demonstrates sweeping a circular profile along a helical path." << endl;
	cout << endl;

	// Create a profile curve: a circle
	NurbsCurve profile;
	NurbsFactory::create_circle(1.0, profile);
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
	NurbsUtil::create_curve_from_points(pathPoints, 3, path);
	cout << "Created helical path with " << numPoints << " control points" << endl;
	cout << endl;

	// Create a closed sweep solid
	cout << "Creating sweep solid with disk caps..." << endl;
	NurbsSolid solid;
	NurbsSweep::sweep_solid(profile, path, solid);

	// Export to OBJ
	cout << "Exporting to OBJ..." << endl;
	OBJWriter ow;
	Mesh mesh;
	NurbsUtil::to_mesh(solid, mesh, 5);
	ow.open("sample_nurbs_sweep.obj");
	ow.write(mesh);
	ow.close();

	// Export both to STEP
	cout << "Exporting to STEP..." << endl;
	StepWriter sw;
	sw.open("sample_nurbs_sweep.step");
	sw.write(solid);
	sw.close();

	cout << "Done!" << endl;
	return 0;
}