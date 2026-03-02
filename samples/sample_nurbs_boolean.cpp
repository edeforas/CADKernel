#include "NurbsBoolean.h"
#include "NurbsBooleanFallback.h"
#include "NurbsIntersection.h"

#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsTrimmedSurface.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "Transform.h"
#include "OBJFile.h"
#include "Mesh.h"
#include "MeshBoolean.h"
#include <iostream>
using namespace std;

int main()
{
	cout << endl << "running sample_nurbs_boolean: compute and save torus-sphere difference" << endl;

	NurbsSolid torus, sphere,diff;
	NurbsFactory::create_torus(8., 2., torus);
	NurbsFactory::create_sphere(8., sphere);
	Translation t(0., 3., 0.);
	sphere.apply_transform(t);
	NurbsBoolean nb;
	nb.compute_difference(torus,sphere,diff);
	
	StepWriter sw;
	sw.open("sample_nurbsboolean.step");
	sw.write(diff);
	sw.close();

	Mesh diffMesh;
	NurbsUtil::to_mesh(diff,diffMesh);
	OBJWriter o;
	o.open("sample_nurbsboolean.obj");
	o.write(diffMesh);
	o.close();

	Mesh torusMesh, sphereMesh;
	NurbsUtil::to_mesh(torus, torusMesh);
	NurbsUtil::to_mesh(sphere, sphereMesh);
	MeshBoolean mb;
	Mesh diffMeshApprox;
	mb.compute_difference(torusMesh, sphereMesh, diffMeshApprox);
	o.open("sample_meshboolean.obj");
	o.write(diffMeshApprox);
	o.close();

	NurbsBooleanFallback fallback;
	std::vector<NurbsTrimmedSurface> trimmedResult;
	if (nb.boolean_difference_trimmed_exact(torus, sphere, trimmedResult))
	{
		Mesh diffMeshFallback;
		NurbsUtil::to_mesh(trimmedResult, diffMeshFallback);
		o.open("sample_nurbsboolean_fallback.obj");
		o.write(diffMeshFallback);
		o.close(); 
	} 
	else
	{
		Mesh diffMeshFallback;
		fallback.difference_mesh(torus, sphere, diffMeshFallback);
		o.open("sample_nurbsboolean_fallback.obj");
		o.write(diffMeshFallback);
		o.close();
	}

	return 0;
}
