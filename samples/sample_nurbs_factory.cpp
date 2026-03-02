
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "Transform.h"
#include "OBJFile.h"

#include <iostream>
using namespace std;

int main()
{
	StepWriter sw;
	sw.open("sample_nurbs_factory.step");

	OBJWriter o;
	Mesh m;
	o.open("sample_nurbs_factory.obj");

	cout << "Creating torus" << endl;
	NurbsSolid torus;
	NurbsFactory::create_torus(10., 3., torus);
	NurbsUtil::to_mesh(torus, m, 8);
	sw.write(torus);
	o.write(m);

	cout << "Creating box" << endl;
	NurbsSolid box;
	NurbsFactory::create_box(10., 8., 6., box);
	box.apply_transform(Translation(Point3(20., 0., 0.)));
	NurbsUtil::to_mesh(box, m, 8);
	sw.write(box);
	o.write(m);

	cout << "Creating sphere" << endl;
	NurbsSolid sphere;
	NurbsFactory::create_sphere(10., sphere);
	sphere.apply_transform(Translation(Point3(40., 0., 0.)));
	NurbsUtil::to_mesh(sphere, m, 8);
	sw.write(sphere);
	o.write(m);


	sw.close();
	o.close();
	return 0;
}
