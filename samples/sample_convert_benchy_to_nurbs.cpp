#include "Mesh.h"

#include "MeshUtil.h"
#include "MeshSolid.h"
#include "MeshSolidUtil.h"
#include "STLFile.h"
#include "OBJFile.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "NurbsSolid.h"
#include <iostream>
using namespace std;

int main()
{
    // Load the Benchy mesh from an OBJ file
    Mesh benchyMesh;
    STLFile::load("3dbenchy.stl", benchyMesh);

    cout << "convert mesh to MeshSolid" << endl;
    MeshSolid benchyMeshSolid;
    MeshSolidUtil::generate_faces(benchyMesh, benchyMeshSolid);
    OBJFile::save("3dbenchy.obj", benchyMeshSolid);

    cout << "Convert the MeshSolid to a NURBS solid" << std::endl;
    NurbsSolid ns;
    NurbsUtil::convert_mesh_to_nurbs(benchyMeshSolid, ns);

    cout << "Saving the NURBS solid to a .step file" << endl;
    StepWriter sw;
    sw.open("3dbenchy.step");
    sw.write(ns);
    sw.close();

    return 0;   
}