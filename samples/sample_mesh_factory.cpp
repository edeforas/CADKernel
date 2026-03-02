#include "MeshFactory.h"
#include "Mesh.h"
#include "OBJFile.h"
#include "Transform.h"

///////////////////////////////////////////////////////////////////////////
int main()
{
	OBJWriter ow;
	ow.open("sample_meshfactory.obj");

	Mesh m;
	MeshFactory::create_tetrahedron(25.,m);
	ow.write(m);

	MeshFactory::create_box(40., 50., 60., m);
	m.apply_transform(Translation(Point3(100., 0., 0.)));
	ow.write(m);

	MeshFactory::create_cylinder(60., 30.,10, m);
	m.apply_transform(Translation(Point3(200., 0., 0.)));
	ow.write(m);

	MeshFactory::create_octahedron(40.,m);
	m.apply_transform(Translation(Point3(0., 100., 0.)));
	ow.write(m);

	MeshFactory::create_dodecahedron(20.,m);
	m.apply_transform(Translation(Point3(100., 100., 0.)));
	ow.write(m);

	MeshFactory::create_icosahedron(30.,m);
	m.apply_transform(Translation(Point3(200, 100, 0.)));
	ow.write(m);

	MeshFactory::create_sphere_uv(30.,8,m);
	m.apply_transform(Translation(Point3(0., 200., 0)));
	ow.write(m);

	MeshFactory::create_sphere_geodesic(30.,10,m);
	m.apply_transform(Translation(Point3(100., 200., 0)));
	ow.write(m);

	MeshFactory::create_torus(25.,15.,8,m);
	m.apply_transform(Translation(Point3(200., 200., 0.)));
	ow.write(m);

	ow.close();
	return 0;
}
///////////////////////////////////////////////////////////////////////////