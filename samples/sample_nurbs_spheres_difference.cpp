#include "NurbsBoolean.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsUtil.h"
#include "OBJFile.h"
#include "StepFile.h"
#include "Transform.h"

#include <iostream>
#include <vector>

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
	std::vector<NurbsTrimmedSurface> exactTrimmedDifference;
	const bool ok = nb.compute_difference(sphereA, sphereB, exactTrimmedDifference, nullptr);
	if (!ok || exactTrimmedDifference.empty())
	{
		std::cerr << "Exact NURBS difference failed for overlapping spheres." << std::endl;
		return 1;
	}

	StepWriter swA;
	swA.open("sample_nurbs_spheres_A.step");
	swA.write(sphereA);
	swA.close();

	StepWriter swB;
	swB.open("sample_nurbs_spheres_B.step");
swB.write(sphereB);
	swB.close();

	StepWriter sw;
	sw.open("sample_nurbs_spheres_difference.step");
	for (const auto& ts : exactTrimmedDifference)
		sw.write(ts);
	sw.close();

	Mesh resultMesh;
	NurbsUtil::to_mesh(exactTrimmedDifference, resultMesh, 24);
	Mesh sphereAMesh;
	Mesh sphereBMesh;
	NurbsUtil::to_mesh(sphereA, sphereAMesh, 24);
	NurbsUtil::to_mesh(sphereB, sphereBMesh, 24);

	OBJWriter ow;
	ow.open("sample_nurbs_spheres_difference.obj");
	ow.write(resultMesh);
	ow.close();

	OBJWriter owA;
	owA.open("sample_nurbs_spheres_A.obj");
	owA.write(sphereAMesh);
	owA.close();

	OBJWriter owB;
	owB.open("sample_nurbs_spheres_B.obj");
	owB.write(sphereBMesh);
	owB.close();

	// Offset instances to compare inputs and exact result side-by-side in one OBJ.
	sphereAMesh.apply_transform(Translation(-28., 0., 0.));
	sphereBMesh.apply_transform(Translation(0., 0., 0.));
	resultMesh.apply_transform(Translation(28., 0., 0.));

	Mesh allMeshes;
	allMeshes.add_mesh(sphereAMesh);
	allMeshes.add_mesh(sphereBMesh);
	allMeshes.add_mesh(resultMesh);
	OBJWriter owAll;
	owAll.open("sample_nurbs_spheres_all.obj");
	owAll.write(allMeshes);
	owAll.close();

	return 0;
}