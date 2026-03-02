
#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "StepFile.h"
#include "Transform.h"
#include "OBJFile.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

// return a random double between min and max, but the same for the same point p
// so that the bump is consistent across seams
double spatial_random_double(const Point3& p,double min, double max) 
{
	double d=cos(p.x()*1372761.12+p.y()*175815.7128+p.z()*265273.76+97.44);
	return min + (max - min) * (d+1.)/2.;
}

int main()
{
	cout << endl << "running sample_bump_cube: create a cube, bump faces, and save" << endl;

	NurbsSolid cube;
	NurbsFactory::create_box(10., 10., 10., cube);

	for (auto& s : cube.surfaces())
	{
		//set to degree 3
		s.degree_elevation_u();
		s.degree_elevation_u();
		s.degree_elevation_v();
		s.degree_elevation_v();

		// add more control points
		for (int i = 1; i < 7; i++)
			s.insert_knot_u(i / 7.);

		for (int j = 1; j < 7; j++)
			s.insert_knot_v(j / 7.);

		//perturbate control points along center direction to create a bump
		for (int i = 0; i < s.nb_points() ; i++)
			{
				auto& p = s.points()[i];
				auto direction = p.normalized();
				p += direction * spatial_random_double(p,-1., 1.);
			}
	}

	StepWriter sw;
	sw.open("sample_bump_cube.step");
	sw.write(cube);
	sw.close();

	OBJWriter o;
	Mesh m;
	NurbsUtil::to_mesh(cube, m, 8);
	o.open("sample_bump_cube.obj");
	o.write(m);
	o.close();

	return 0;
}
