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
	cout << "This sample demonstrates sweeping a circular profile along a sinusoidal path." << endl;
	cout << "It creates both a standard sweep and a perpendicular sweep for comparison." << endl;
	cout << endl;

	// Create a profile curve: a circle
	NurbsCurve profile;
	NurbsFactory::create_circle(2.0, profile); // radius 2
	cout << "Created circular profile with radius 2.0" << endl;

	// Create a path curve: a smooth sinusoidal path
	NurbsCurve path;
	std::vector<Point3> pathPoints;
	int numPoints = 20;
	for (int i = 0; i < numPoints; ++i) {
		double t = (double)i / (numPoints - 1);
		double x = 0.0;
		double y = 5.0 * sin(t * 2.0 * acos(-1.0)); // Sinusoidal in y
		double z = t * 10.0;
		pathPoints.push_back(Point3(x, y, z));
	}
	NurbsUtil::create_curve_from_points(pathPoints, 3, path); // degree 3 for smooth curve
	cout << "Created sinusoidal path with " << numPoints << " control points" << endl;

	// Sweep the profile along the path to create a surface (standard sweep)
	cout << "Performing standard sweep (profile orientation fixed)..." << endl;
	NurbsSurface surface;
	bool success = NurbsSweep::sweep(profile, path, surface);
	if (!success) {
		cout << "Failed to sweep surface!" << endl;
		return 1;
	}

	// Also create a perpendicular sweep
	cout << "Performing perpendicular sweep (profile stays perpendicular to path)..." << endl;
	NurbsSurface surfacePerp;
	bool successPerp = NurbsSweep::sweep(profile, path, surfacePerp, true);
	if (!successPerp) {
		cout << "Failed to sweep perpendicular surface!" << endl;
		return 1;
	}

	// Force non-closed to avoid STEP issues
	//surface.set_closed_u(false);
	//surface.set_closed_v(false);
	//surfacePerp.set_closed_u(false);
	//surfacePerp.set_closed_v(false);

	// Save to STEP files
	cout << "Saving to STEP files..." << endl;
	StepWriter sw;
	sw.open("sample_nurbs_sweep.step");
	sw.write(surface);
	sw.close();

	sw.open("sample_nurbs_sweep_perp.step");
	sw.write(surfacePerp);
	sw.close();

	// Save to OBJ files
	cout << "Saving to OBJ files..." << endl;
	Mesh mesh;
	NurbsUtil::to_mesh(surface, mesh, 5);
	OBJWriter ow;
	ow.open("sample_nurbs_sweep.obj");
	ow.write(mesh);
	ow.close();

	Mesh meshPerp;
	NurbsUtil::to_mesh(surfacePerp, meshPerp, 5);
	ow.open("sample_nurbs_sweep_perp.obj");
	ow.write(meshPerp);
	ow.close();

	cout << endl;
	cout << "Sample completed successfully!" << endl;
	cout << "Files created:" << endl;
	cout << "  - sample_nurbs_sweep.step      (standard sweep)" << endl;
	cout << "  - sample_nurbs_sweep.obj       (standard sweep)" << endl;
	cout << "  - sample_nurbs_sweep_perp.step (perpendicular sweep)" << endl;
	cout << "  - sample_nurbs_sweep_perp.obj  (perpendicular sweep)" << endl;
	cout << endl;
	cout << "The perpendicular sweep keeps the profile circle perpendicular to the path" << endl;
	cout << "direction at each point, creating a more natural tube-like surface." << endl;

	return 0;
}