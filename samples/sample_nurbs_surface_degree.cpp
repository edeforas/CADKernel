#include "NurbsSurface.h"
#include "NurbsUtil.h"

#include "OBJFile.h"
#include "StepFile.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

///////////////////////////////////////////////////////////////////////////
int main()
{
	cout << endl << "sample_nurbs_surface_degree" << endl;
	cout << "Create several nurbsurfaces with degree from 1 to 3 in X and Y" << endl;
	
	int nbPointsU = 9;
	int nbPointsV = 7;

	// create surface
	vector<Point3> points,pt;
	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
			points.push_back(Point3(u, v, 2.*(float)rand() / RAND_MAX));

	OBJWriter ow;
	ow.open("sample_nurbs_surface_degree.obj");

	StepWriter sw;
	sw.open("sample_nurbs_surface_degree.step");

	for(int uDeg=1; uDeg<=3;uDeg++)
		for (int vDeg = 1; vDeg <= 3; vDeg++)
		{
			NurbsSurface n;

			n.set_degree(uDeg, vDeg);

			//translate point
			pt = points;
			for (auto& p : pt)
			{
				p.x() += uDeg * nbPointsU*1.1;
				p.y() += vDeg * nbPointsV*1.1;
			}

			n.set_points(pt, nbPointsU, nbPointsV);
			n.set_uniform_u();
			n.set_uniform_v();
			n.set_equals_weights();

			Mesh m;
			NurbsUtil::to_mesh(n, m, 10);

			ow.write(m);
			sw.write(n);
		}
	ow.close();
	return 0;
}
