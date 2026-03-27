
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
	NurbsFactory::create_torus(6., 3., torus);
	NurbsUtil::to_mesh(torus, m, 8);
	sw.write(torus);
	o.write(m);

	cout << "Creating box" << endl;
	NurbsSolid box;
	NurbsFactory::create_box(11., 9., 10., box);
	box.apply_transform(Translation(Point3(20., 0., 0.)));
	NurbsUtil::to_mesh(box, m, 8);
	sw.write(box);
	o.write(m);

	cout << "Creating sphere" << endl;
	NurbsSolid sphere;
	NurbsFactory::create_sphere(8, sphere);
	sphere.apply_transform(Translation(Point3(40., 0., 0.)));
	NurbsUtil::to_mesh(sphere, m, 8);
	sw.write(sphere);
	o.write(m);

	cout << "Creating cylinder" << endl;
	NurbsSolid cylinder;
	NurbsFactory::create_cylinder(7, 10, cylinder);
	cylinder.apply_transform(Translation(Point3(0., 20., 0.)));
	NurbsUtil::to_mesh(cylinder, m, 8);
	sw.write(cylinder);
	o.write(m);

	sw.close();
	o.close();
	return 0;
}
