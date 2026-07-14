#include "Mesh.h"

#include "MeshUtil.h"
#include "MeshSolid.h"
#include "MeshSolidUtil.h"
#include "STLFile.h"
#include "OBJFile.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "NurbsSolid.h"

int main()
{
    // Load the Benchy mesh from an OBJ file
    Mesh benchyMesh;
    STLFile::load("3dbenchy.stl", benchyMesh);

    // convert mesh to MeshSolid
    MeshSolid benchyMeshSolid;
    MeshSolidUtil::generate_faces(benchyMesh, benchyMeshSolid);
    OBJFile::save("3dbenchy.obj", benchyMeshSolid);

    // Convert the MeshSolid to a NURBS solid
    NurbsSolid ns;
    NurbsUtil::convert_mesh_to_nurbs(benchyMeshSolid, ns);

    // Save the NURBS solid to a file
    StepWriter sw;
    sw.open("3dbenchy.step");
    sw.write(ns);
    sw.close();

    return 0;   
}