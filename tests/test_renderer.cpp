#include "Renderer.h"
#include "Image.h"
#include "ImageIoBmp.h"
#include "MeshFactory.h"
#include "NurbsCurve.h"
#include "NurbsFactory.h"
#include "NurbsUtil.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "Transform.h"

#include <iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
int main()
{
	int iWidth = 1600;
	int iHeight = 1200;
	double dAngleX = 20., dAngleY = 10., dAhead = 50., dZoom = 2000.;
	
	Mesh mTorus;
	MeshFactory::create_torus(15, 3, 32,mTorus);
	mTorus.set_color(GREY);

	Mesh mSphere;
	MeshFactory::create_sphere_geodesic(5.,6, mSphere);
	mSphere.set_color(DARK_GREEN);
	mSphere.apply_transform(Translation(Point3(10, 0., 0.)));

	NurbsSurface nSurface;
	NurbsFactory::create_torus(8., 2., nSurface);
	for (auto& p : nSurface.points())
		p += Point3(-12., 0., 0.);

	NurbsSolid nSolid;
	NurbsFactory::create_box(10., 10., 10., nSolid);
	nSolid.apply_transform(Translation(Point3(0., -12., 0.)));

	Image img(iWidth, iHeight, 4);
	Renderer eng((int*)img.data(),iWidth, iHeight);
	eng.set_background(DARK_GREY);
	eng.add_ambient_light(GREY, 1.);
	eng.add_diffuse_light(RED, 1., Point3(1., 0., 1.));
	eng.add_diffuse_light(GREEN, 1., Point3(0., 1., 0.));
	eng.add_diffuse_light(BLUE, 1., Point3(0., 1., 1.));

	NurbsCurve n;
	int degree = 2;
	vector<Point3> points,nurbPL;
	int nbPoint = 5;

	for (int i = 0; i < nbPoint; i++)
	{
		Point3 p(Point3((double)rand() / RAND_MAX-0.5, (double)rand() / RAND_MAX-0.5, (double)rand() / RAND_MAX-0.5));
		points.push_back(p*30.);
	}

	NurbsUtil::create_curve_from_points(points, degree,n);
	n.to_polyline(nurbPL);

	for (int i = 0; i < 360; i += 40)
	{
		Mesh m;
		eng.clear();
		eng.set_camera(0., 0., 0., dAhead, i, i/2., i/3., dZoom);
		
		eng.draw_mesh(mTorus,false);
		eng.draw_mesh(mSphere,true);
		eng.draw_surface(nSurface, 5, false, BLUE);
		eng.draw_solid(nSolid, 4, true, GREY);
		eng.draw_polyline(nurbPL, WHITE);
		eng.draw_polyline(points, GREEN);

		string sFile = string("rendered_") + to_string(i) + ".bmp";

		cout << "Writing: " << sFile << endl;
		if (!ImageIoBmp::write(sFile, &img))
		{
			cout << "Unable to save rendered image" << sFile << endl;
			return -1;
		}
	}

	cout << "Test Finished.";
	return 0;
}
/////////////////////////////////////////////////////////////////////////////