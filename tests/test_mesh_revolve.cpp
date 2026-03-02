#include "MeshFactory.h"
#include "Geometry.h"

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

void test_mesh_revolve_basic()
{
	std::vector<Point3> profile;
	profile.push_back(Point3(10., 0., -5.));
	profile.push_back(Point3(10., 0., 0.));
	profile.push_back(Point3(10., 0., 5.));

	Mesh m;
	const int nbSegments = 12;
	MeshFactory::create_revolve(profile, nbSegments, m);

	test_bool(m.nb_vertices() == nbSegments * (int)profile.size(), "unexpected revolved mesh vertex count");
	test_bool(m.nb_triangles() == nbSegments * ((int)profile.size() - 1) * 2, "unexpected revolved mesh triangle count");

	for (int i = 0; i < m.nb_vertices(); ++i)
	{
		Point3 p;
		m.get_vertex(i, p);
		double radius = std::sqrt(p.x() * p.x() + p.y() * p.y());
		test_bool(std::fabs(radius - 10.) < 1.e-9, "revolved vertex has unexpected radius");
	}
}

void test_mesh_revolve_capped()
{
	std::vector<Point3> profile;
	profile.push_back(Point3(8., 0., -4.));
	profile.push_back(Point3(8., 0., 0.));
	profile.push_back(Point3(8., 0., 4.));

	Mesh m;
	const int nbSegments = 10;
	MeshFactory::create_revolve(profile, nbSegments, true, m);

	const int expectedVertices = nbSegments * (int)profile.size() + 2;
	const int expectedTriangles = nbSegments * ((int)profile.size() - 1) * 2 + 2 * nbSegments;
	test_bool(m.nb_vertices() == expectedVertices, "unexpected capped revolved mesh vertex count");
	test_bool(m.nb_triangles() == expectedTriangles, "unexpected capped revolved mesh triangle count");
}

int main()
{
	test_mesh_revolve_basic();
	test_mesh_revolve_capped();
	cout << "Test Finished.";
	return 0;
}
