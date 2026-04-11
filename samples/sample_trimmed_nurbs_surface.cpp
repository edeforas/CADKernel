#include "NurbsFactory.h"
#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsUtil.h"
#include "Mesh.h"
#include "OBJFile.h"
#include "StepFile.h"
#include "Transform.h"
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

///////////////////////////////////////////////////////////////////////////
int main()
{
	cout << endl << "sample_trimmed_nurbs_surface" << endl;

	// Create a simple rectangular NURBS surface
	NurbsSurface base;
	NurbsFactory::create_quad(
		Point3(-2., -2., 0.),    // bottom-left
		Point3(2., -2., 0.),   // bottom-right
		Point3(2., 2., 0.),  // top-right
		Point3(-2., 2., 0.),   // top-left
		base);

	// Create a trimmed surface from the base
	NurbsTrimmedSurface trimmedSurface(base);

	// Add outer loop (boundary of the surface in UV space)
	trimmedSurface.add_full_outer_loop();

	// Create a circular hole in UV parameter space
	// The surface domain is [0,1] x [0,1] in UV coordinates
	// Create circle with radius 0.15 centered at UV (0.5, 0.5)
	NurbsCurve circularHole;
	NurbsFactory::create_circle(0.15, circularHole);
	
	// The create_circle creates a 3D circle. For a 2D UV curve, we keep X,Y as U,V
	// and make sure it's in [0,1]^2 parameter space
	// Translate center to (0.5, 0.5) in UV space
	Translation(Point3(0.5, 0.5, 0.)).apply_all(circularHole.points());
	trimmedSurface.add_inner_loop(circularHole);

	// Convert to mesh for visualization
	Mesh mesh;
	NurbsUtil::to_mesh(trimmedSurface, mesh, 40);
	
	// Save as OBJ file
	OBJWriter objWriter;
	objWriter.open("sample_trimmed_nurbs_surface.obj");
	objWriter.write(mesh);
	objWriter.close();
	cout << "Saved: sample_trimmed_nurbs_surface.obj" << endl;

	// Save as STEP file
	StepWriter stepWriter;
	stepWriter.open("sample_trimmed_nurbs_surface.step");
	stepWriter.write(trimmedSurface);
	stepWriter.close();
	cout << "Saved: sample_trimmed_nurbs_surface.step" << endl;

	cout << "Sample completed successfully!" << endl << endl;

	return 0;
}
