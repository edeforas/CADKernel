#include "NurbsBooleanFallback.h"

#include "MeshBoolean.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"

namespace
{
	void solids_to_meshes(const NurbsSolid& a, const NurbsSolid& b, Mesh& meshA, Mesh& meshB, int iNbSegments)
	{
		if (iNbSegments < 4)
			iNbSegments = 4;

		NurbsUtil::to_mesh(a, meshA, iNbSegments);
		NurbsUtil::to_mesh(b, meshB, iNbSegments);
	}
}

void NurbsBooleanFallback::union_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments)
{
	Mesh meshA, meshB;
	solids_to_meshes(a, b, meshA, meshB, iNbSegments);

	MeshBoolean mb;
	mb.compute_union(meshA, meshB, result);
}

void NurbsBooleanFallback::intersection_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments)
{
	Mesh meshA, meshB;
	solids_to_meshes(a, b, meshA, meshB, iNbSegments);

	MeshBoolean mb;
	mb.compute_intersection(meshA, meshB, result);
}

void NurbsBooleanFallback::difference_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments)
{
	Mesh meshA, meshB;
	solids_to_meshes(a, b, meshA, meshB, iNbSegments);

	MeshBoolean mb;
	mb.compute_difference(meshA, meshB, result);
}