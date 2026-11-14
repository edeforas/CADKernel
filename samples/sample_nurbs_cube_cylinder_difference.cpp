#include "NurbsBoolean.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "OBJFile.h"
#include <iostream>
using namespace std;

int main()
{
	cout << endl << "running sample_nurbs_cube_cylinder_difference: cube minus cylinder" << endl;

	NurbsSolid box, cylinder;
	NurbsFactory::create_box(20.0, 20.0, 20.0, box);
	NurbsFactory::create_cylinder(6.0, 30.0, cylinder);

	// Align cylinder through the center of the box
	// Both shapes are centered at origin by default

	NurbsBoolean nb;
	NurbsSolid difference;
	nb.compute_difference(box, cylinder, difference);
	
	// Save resulting solid to STEP
	StepWriter sw;
	sw.open("sample_nurbs_cube_cylinder_difference.step");
	sw.write(difference);
	sw.close();

	// Convert result to mesh and save to OBJ
	Mesh resultMesh;
	NurbsUtil::to_mesh(difference, resultMesh, 32);

	OBJWriter ow;
	ow.open("sample_nurbs_cube_cylinder_difference.obj");
	ow.write(resultMesh);
	ow.close();

	cout << "Wrote sample_nurbs_cube_cylinder_difference.step and sample_nurbs_cube_cylinder_difference.obj" << endl;
	return 0;
}
