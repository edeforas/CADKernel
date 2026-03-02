#include "MeshFactory.h"

#include "STLFile.h"
#include "OBJFile.h"

///////////////////////////////////////////////////////////////////////////
int main()
{
	// build object
	Mesh m;
	MeshFactory::create_cylinder(100., 40., 8, m);

	//test save
	STLFile::save("cylinder.stl", m);
	OBJFile::save("cylinder.obj", m);

	//test load and save stl
	Mesh m1;
	STLFile::load("cylinder.stl", m1);
	STLFile::save("cylinder.stl", m1);

	//test load and save obj
	Mesh m2;
	OBJFile::load("cylinder.obj", m2);
	OBJFile::save("cylinder.obj", m2);

	return 0;
}
///////////////////////////////////////////////////////////////////////////