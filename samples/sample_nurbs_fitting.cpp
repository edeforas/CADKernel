#include "NurbsFitting.h"
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
	cout << "NURBS Fitting Sample" << endl;
	cout << "Creating sample points for a wavy surface..." << endl;

	// Create a grid of sample points for a wavy surface
	const int numSamplesU = 8;
	const int numSamplesV = 8;
	const double uMin = -5.0, uMax = 15.0;
	const double vMin = -5.0, vMax = 5.0;

	std::vector<Point3> samplePoints;
	
	for (int j = 0; j < numSamplesV; ++j) {
		double v = vMin + (vMax - vMin) * j / (numSamplesV - 1);
		for (int i = 0; i < numSamplesU; ++i) {
			double u = uMin + (uMax - uMin) * i / (numSamplesU - 1);

			// Create a wavy surface: z = sin(u) * cos(v) + 0.5 * sin(2*u)
			double z = sin(u) * cos(v) + 0.5 * sin(2.0 * u)+u*u/10.;

			samplePoints.push_back(Point3(u, v, z));
		}
	}

	cout << "Generated " << samplePoints.size() << " sample points" << endl;

	// Fit a NURBS surface to the sample points
	cout << "Fitting NURBS surface..." << endl;
	NurbsSurface fittedSurface;

	const int degreeU = 3;
	const int degreeV = 3;
	const int numCtrlU = 8;
	const int numCtrlV = 8;

	bool success = NurbsFitting::fit_surface_least_squares(
		samplePoints,
		numSamplesU,
		numSamplesV,
		degreeU,
		degreeV,
		numCtrlU,
		numCtrlV,
		fittedSurface
	);

	if (!success) {
		cout << "Failed to fit NURBS surface!" << endl;
		return 1;
	}

	cout << "Successfully fitted NURBS surface with " << fittedSurface.nb_points_u()
		 << "x" << fittedSurface.nb_points_v() << " control points" << endl;

	// Save to STEP file
	cout << "Saving to STEP file..." << endl;
	StepWriter sw;
	sw.open("sample_nurbs_fitting.step");
	sw.write(fittedSurface);
	sw.close();

	// Save to OBJ file
	cout << "Saving to OBJ file..." << endl;
	Mesh mesh;
	NurbsUtil::to_mesh(fittedSurface, mesh, 32); // High resolution for smooth surface

	OBJWriter ow;
	ow.open("sample_nurbs_fitting.obj");
	ow.write(mesh);
	ow.close();

	cout << "Sample completed successfully!" << endl;
	cout << "Files created:" << endl;
	cout << "  - sample_nurbs_fitting.step" << endl;
	cout << "  - sample_nurbs_fitting.obj" << endl;

	return 0;
}