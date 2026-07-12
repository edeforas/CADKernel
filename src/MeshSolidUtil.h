#ifndef MeshSolidUtil_
#define MeshSolidUtil_

#include <vector>
using namespace std;

#include "Geometry.h"
#include "MeshSolid.h"

class Mesh;
///////////////////////////////////////////////////////////////////////////
namespace MeshSolidUtil
{
	void generate_faces(const Mesh& m,  MeshSolid& ms, double dToleranceAngleDeg=0.1);

}
///////////////////////////////////////////////////////////////////////////

#endif