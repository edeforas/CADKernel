#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"

#include "StepFile.h"
#include "STLFile.h"

#include <iostream> 
using namespace std;


int main()
{
	cout << endl << "simple sample_cube: create a cube, save as .step and .stl file" << endl;

	NurbsSolid cube;
	NurbsFactory::create_box(10., 10., 10., cube);

	StepWriter sw;
	sw.open("sample_cube.step");
	sw.write(cube);
	sw.close();

	Mesh m;
	NurbsUtil::to_mesh(cube, m);
	STLFile::save("sample_cube.stl", m);

	return 0;
}