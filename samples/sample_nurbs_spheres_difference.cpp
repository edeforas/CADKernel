#include "NurbsBoolean.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "OBJFile.h"
#include "StepFile.h"
#include "Transform.h"

#include <iostream>


int main()
{
	std::cout << "running sample_nurbs_spheres_difference" << std::endl;

	NurbsSolid sphereA;
	NurbsSolid sphereB;
	NurbsFactory::create_sphere(10., sphereA);
	NurbsFactory::create_sphere(10., sphereB);

	// Partial overlap so the boolean result has a non-trivial exact trim boundary.
	Translation moveB(7., 0., 0.);
	sphereB.apply_transform(moveB);

	NurbsBoolean nb;
	NurbsSolid solidDifference;
	nb.compute_difference(sphereA, sphereB, solidDifference);

	StepWriter swStep;
	swStep.open("sample_nurbs_spheres_difference.step");
	swStep.write(solidDifference);
	swStep.close();

	Mesh resultMesh;
	NurbsUtil::to_mesh(solidDifference, resultMesh, 24);
	OBJWriter ow;
	ow.open("sample_nurbs_spheres_difference.obj");
	ow.write(resultMesh);
	ow.close();

	return 0;
}