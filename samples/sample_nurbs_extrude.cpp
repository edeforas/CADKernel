#include "NurbsExtrude.h"
#include "NurbsFactory.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "OBJFile.h"
#include "Mesh.h"

#include <iostream>
#include <vector>
using namespace std;

int main()
{
	cout << "NURBS Extrude Sample" << endl;

	// Create a NURBS cube
	NurbsSolid cube;
	NurbsFactory::create_box(10., 10., 10., cube);

	cout << "Original cube has " << cube.surfaces().size() << " surfaces" << endl;

	// Compute extrusion directions using face normals at center (u=0.5, v=0.5)
	vector<Point3> directions;
	const auto& surfaces = cube.surfaces();
	for (const auto& surface : surfaces) {
		Point3 normal;
		if (surface.normal(0.5, 0.5, normal)) {
			// Extrude in the direction of the normal, scaled by 5 units
			directions.push_back(normal * 5.0);
		} else {
			cout << "Failed to compute normal for a surface" << endl;
			directions.push_back(Point3(0, 0, 5)); // fallback
		}
	}

	// Extrude all faces
	NurbsExtrude extruder;
	int extrudedCount = 0;
	size_t originalSurfaceCount = surfaces.size();
	for (size_t i = 0; i < originalSurfaceCount; ++i) {
		bool success = extruder.extrude_face(cube, i, directions[i]);
		if (success) {
			extrudedCount++;
		} else {
			cout << "Failed to extrude face " << i << endl;
		}
	}

	cout << "Successfully extruded " << extrudedCount << " faces" << endl;
	cout << "Cube now has " << cube.surfaces().size() << " surfaces" << endl;

	// Save to STEP file
	cout << "Saving to STEP file..." << endl;
	StepWriter sw;
	sw.open("sample_nurbs_extrude.step");
	sw.write(cube);
	sw.close();

	// Save to OBJ file
	cout << "Saving to OBJ file..." << endl;
	Mesh mesh;
	NurbsUtil::to_mesh(cube, mesh, 16);
	OBJWriter ow;
	ow.open("sample_nurbs_extrude.obj");
	ow.write(mesh);
	ow.close();

	cout << "Sample completed successfully!" << endl;
	cout << "Files created:" << endl;
	cout << "  - sample_nurbs_extrude.step" << endl;
	cout << "  - sample_nurbs_extrude.obj" << endl;

	return 0;
}