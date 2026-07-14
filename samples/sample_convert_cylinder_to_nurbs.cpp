#include "Mesh.h"

#include "MeshUtil.h"
#include "MeshSolid.h"
#include "MeshSolidUtil.h"
#include "MeshFactory.h"
#include "OBJFile.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "NurbsSolid.h"
#include <iostream>
using namespace std;

int main()
{
    // Create a cylinder without face hierarchy
    Mesh m;
    MeshFactory::create_cylinder(30,10,8,m);

    cout << "convert mesh to MeshSolid" << endl;
    MeshSolid mSolid;
    MeshSolidUtil::generate_faces(m, mSolid,12.0);
    OBJFile::save("cylinder.obj", mSolid);

    cout << "Convert the MeshSolid to a Nurbs solid" << std::endl;
    NurbsSolid ns;
    NurbsUtil::convert_mesh_to_nurbs(mSolid, ns);

    cout << "Saving the Nurbs solid to a .step file" << endl;
    StepWriter sw;
    sw.open("cylinder_from_mesh.step");
    sw.write(ns);
    sw.close();

    return 0;   
}