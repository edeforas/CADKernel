#include "STLFile.h"
#include "Mesh.h"

#include <vector>
#include <iostream>
using namespace std;

struct Tetrahedron
{
	Point3 p1, p2, p3, p4;
};

///////////////////////////////////////////////////////////////////////////
int main()
{
	//create Shierpinsky Fractal tetrahedron by subdivision
	cout << "Fractal computation Start." << endl;
	
	//create first plain tetrahedron
	Tetrahedron t;
	t.p1 = Point3( 1.,  1.,  1.);
	t.p2 = Point3(-1., -1.,  1.);
	t.p3 = Point3( 1., -1., -1.);
	t.p4 = Point3(-1.,  1., -1.);
	vector<Tetrahedron> vt;
	vt.push_back(t);

	//subdivide tetrahedron in 4 pieces and loop
	for (int iScale = 0; iScale < 3; iScale++)
	{
		vector<Tetrahedron> vt2;
		for (int i = 0; i < vt.size(); i++)
		{
			const auto& t1 = vt[i];

			Tetrahedron t2;
			t2.p1 = t1.p1;
			t2.p2 = (t1.p1 + t1.p2) / 2.;
			t2.p3 = (t1.p1 + t1.p3) / 2.;
			t2.p4 = (t1.p1 + t1.p4) / 2.;
			vt2.push_back(t2);

			t2.p1 = t1.p2;
			t2.p2 = (t1.p2 + t1.p4) / 2.;
			t2.p3 = (t1.p2 + t1.p3) / 2.;
			t2.p4 = (t1.p2 + t1.p1) / 2.;
			vt2.push_back(t2);

			t2.p1 = t1.p3;
			t2.p2 = (t1.p3 + t1.p1) / 2.;
			t2.p3 = (t1.p3 + t1.p2) / 2.;
			t2.p4 = (t1.p3 + t1.p4) / 2.;
			vt2.push_back(t2);

			t2.p1 = t1.p4;
			t2.p2 = (t1.p4 + t1.p2) / 2.;
			t2.p3 = (t1.p4 + t1.p1) / 2.;
			t2.p4 = (t1.p4 + t1.p3) / 2.;
			vt2.push_back(t2);
		}

		vt = vt2;
	}

	//create mesh using smalls tetrahedrons
	Mesh m;
	for (int i = 0; i < vt.size();i++)
	{
		const auto& t1 = vt[i];

		int iV1=m.add_vertex(t1.p1);
		int iV2=m.add_vertex(t1.p2);
		int iV3=m.add_vertex(t1.p3);
		int iV4=m.add_vertex(t1.p4);

		m.add_triangle(iV1, iV2, iV3);
		m.add_triangle(iV1, iV3, iV4);
		m.add_triangle(iV1, iV4, iV2);
		m.add_triangle(iV2, iV4, iV3);
	}

	STLFile::save("sample_fractal_tetrahedron.stl", m);

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////