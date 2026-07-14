#include "Mesh.h"

#include "STLFile.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "NurbsSolid.h"

int main()
{
    // Load the Benchy mesh from an OBJ file
    Mesh benchyMesh;
    STLFile::load("3dbenchy.stl", benchyMesh);

    // Convert the mesh to a NURBS surface
    NurbsSolid ns;
    NurbsUtil::convert_mesh_to_nurbs(benchyMesh, ns);

    // Save the NURBS surface to a file
    StepWriter sw;
    sw.open("3dbenchy.step");
    sw.write(ns);
    sw.close();

    return 0;   
}