#include <cmath>
#include <cassert>
using namespace std;

#include "Renderer.h"
#include "Mesh.h"
#include "NurbsUtil.h"

#define MAXINT (1<<28)

// used by mapfacet
#define my_swap(a,b,t) { t=a; a=b; b=t; }

////////////////////////////////////////////////////////////////////////////////
void Renderer::set_camera(double ox, double oy, double oz, double ahead, double yaw, double pitch, double roll, double zoom)
{
	_camera.set_origin(ox, oy, oz, ahead);
	_camera.set_angles(yaw, pitch, roll);
	_camera.set_screen(_Xmax, _Ymax, zoom);
}
////////////////////////////////////////////////////////////////////////////////
Renderer::Renderer(int* pBuffer, int xm, int ym)
{
	_Xmax = xm;
	_Ymax = ym;
	_pixelBuffer = pBuffer;

	_wBuffer = new float[_Xmax * _Ymax];

	_iBackgroundColor = 0;
	clear();
}
////////////////////////////////////////////////////////////////////////////////
Renderer::~Renderer()
{
	delete[] _wBuffer;

	for (size_t i = 0; i < _lights.size(); i++)
		delete _lights[i];

	_lights.clear();
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::add_ambient_light(int iAmbiantColor, double dAmbiantFactor)
{
	_lights.push_back(new RendererLightAmbiant(iAmbiantColor, dAmbiantFactor));
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::add_diffuse_light(int iDiffuseColor, double dDiffuseFactor, const Point3& direction)
{
	_lights.push_back(new RendererLightDiffuse(iDiffuseColor, dDiffuseFactor, direction));
}
////////////////////////////////////////////////////////////////////////////////
int Renderer::compute_color_with_lights(int iColor, const Point3& normal)
{
	double dRed = ((iColor >> 16) & 0xff) / 256.;
	double dGreen = ((iColor >> 8) & 0xff) / 256.;
	double dBlue = ((iColor) & 0xff) / 256.;

	double dRTotal = 0., dGTotal = 0., dBTotal = 0.;
	for (size_t i = 0; i < _lights.size(); i++)
	{
		double dR = dRed, dG = dGreen, dB = dBlue;
		_lights[i]->apply(dR, dG, dB, normal);
		dRTotal += dR;
		dGTotal += dG;
		dBTotal += dB;
	}

	int iRed = (int)(dRTotal * 256. + 0.5);
	if (iRed > 255)
		iRed = 255;

	int iGreen = (int)(dGTotal * 256. + 0.5);
	if (iGreen > 255)
		iGreen = 255;

	int iBlue = (int)(dBTotal * 256. + 0.5);
	if (iBlue > 255)
		iBlue = 255;

	return (iRed << 16) + (iGreen << 8) + iBlue;
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_mesh(const Mesh& m, bool bDrawEdges, int iColor)
{
	// us default mesh color if iColor is -1
	int color = iColor;
	if (iColor == -1)
		color = m.get_color();

	for (int i = 0; i < m.nb_triangles(); i++)
	{
		if (m.is_triangle_unlinked(i))
			continue;

		Triangle3 t;
		m.get_triangle(i, t);
		bool bDrawed = draw_triangle_1color(t.p1(), t.p2(), t.p3(), color);

		if (bDrawEdges && bDrawed) //todo custom color, clean zbuffer , etc
		{
			draw_line(t.p1(), t.p2(), (0)); //RGB(0,0,0)
			draw_line(t.p1(), t.p3(), (0)); //RGB(0,0,0)
			draw_line(t.p2(), t.p3(), (0)); //RGB(0,0,0)
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_surface(const NurbsSurface& n, int iNbSegments, bool bDrawEdges, int iColor)
{
	Mesh m;
	NurbsUtil::to_mesh(n, m, iNbSegments);
	draw_mesh(m, bDrawEdges, iColor);
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_solid(const NurbsSolid& n, int iNbSegments, bool bDrawEdges, int iColor)
{
	Mesh m;
	NurbsUtil::to_mesh(n, m, iNbSegments);
	draw_mesh(m, bDrawEdges, iColor);
}
////////////////////////////////////////////////////////////////////////////////
bool Renderer::draw_triangle_1color(const Point3& A, const Point3& B, const Point3& C, int iColor, bool bTwofaces)
{
	int ax, ay, bx, by, cx, cy;
	double aw, bw, cw;
	_camera.project(A, ax, ay, aw);
	_camera.project(B, bx, by, bw);
	_camera.project(C, cx, cy, cw);

	// test if seeing the back of the facet
	if (bTwofaces == false)
	{
		Point3 norm2d = (Point3(ax - bx, ay - by, 0.)).cross_product(Point3(ax - cx, ay - cy, 0.));
		if (norm2d.z() <= 0.)
			return false;
	}

	//exit if triangle is trivially out of screen
	if ((ax <= 0.) && (bx <= 0.) && (cx <= 0.))
		return false;
	if ((ay <= 0.) && (by <= 0.) && (cy <= 0.))
		return false;
	if ((ax >= _Xmax) && (bx >= _Xmax) && (cx >= _Xmax))
		return false;
	if ((ay >= _Ymax) && (by >= _Ymax) && (cy >= _Ymax))
		return false;

	//order points in y axis
	if ((by <= ay) && (by <= cy))
	{
		//swap B and A
		int tmpx = ax; ax = bx; bx = tmpx;
		int tmpy = ay; ay = by; by = tmpy;
		double tmpw = aw; aw = bw; bw = tmpw;
	}
	else
		if ((cy <= ay) && (cy <= by))
		{
			//swap C and A
			int tmpx = ax; ax = cx; cx = tmpx;
			int tmpy = ay; ay = cy; cy = tmpy;
			double tmpw = aw; aw = cw; cw = tmpw;
		}

	if (cy < by)
	{
		//swap C and B
		int tmpx = bx; bx = cx; cx = tmpx;
		int tmpy = by; by = cy; cy = tmpy;
		double tmpw = bw; bw = cw; cw = tmpw;
	}

	//split in two trapeze, compute D= the horizontal intersection of B with AC
	double t;
	if (cy != ay)
		t = (double)(by - ay) / (cy - ay);
	else
		t = 0.;

	assert(t >= 0.);
	assert(t <= 1.);

	int dx = (int)(cx * t + (1. - t) * ax + 0.5);
	int dy = (int)(cy * t + (1. - t) * ay + 0.5);
	double dw = cw * t + (1. - t) * aw;

	//reorder in x such that bx<dx
	if ((dx < bx))
	{
		//exchange D and B
		int tmpx = bx; bx = dx; dx = tmpx;
		int tmpy = by; by = dy; //dy = tmpy;
		double tmpw = bw; bw = dw; dw = tmpw;
	}

	assert(bx <= dx);
	assert(ay <= by);
	assert(by <= cy);

	//draw each trapez
	Point3 normal = Triangle3(A, B, C).normal();
	int iColorWithLights = compute_color_with_lights(iColor, normal);
	bool b1 = draw_trapeze(ax, aw, ax, aw, ay, bx, bw, dx, dw, by, iColorWithLights);
	bool b2 = draw_trapeze(bx, bw, dx, dw, by, cx, cw, cx, cw, cy, iColorWithLights);

	return b1 || b2;
}
////////////////////////////////////////////////////////////////////////////////
bool Renderer::draw_trapeze(double ax, double aw, double bx, double bw, int ay, double cx, double cw, double dx, double dw, int cy, int color)
{
	assert(ax <= bx);
	assert(cx <= dx);
	assert(ay <= cy);


	if (cy < 0.)
		return false;

	if (ay >= _Ymax)
		return false;

	//reduce facet y
	if (ay < 0)
	{
		double t = (double)(0. - ay) / (cy - ay);
		assert(t >= 0.);
		assert(t <= 1.);

		ay = 0;
		ax = cx * t + (1. - t) * ax;
		aw = cw * t + (1. - t) * aw;
		bx = dx * t + (1. - t) * bx;
		bw = dw * t + (1. - t) * bw;
	}
	if (cy >= _Ymax)
	{
		double t = (double)(cy - (_Ymax - 1)) / (cy - ay);
		assert(t >= 0.);
		assert(t <= 1.);

		cy = _Ymax - 1;
		cx = ax * t + (1. - t) * cx;
		cw = aw * t + (1. - t) * cw;
		dx = bx * t + (1. - t) * dx;
		dw = bw * t + (1. - t) * dw;
	}

	assert(ay >= 0);
	assert(ay <= cy);
	assert(cy < _Ymax);

	//plot trapeze
	assert(ay <= cy);
	for (int y = ay; y < cy; y++)
	{
		//plot line
		double t = (double)(y - ay) / (cy - ay);
		assert(t >= 0.);
		assert(t <= 1.);

		int p1x = (int)(cx * t + ax * (1. - t) + 0.5);
		double p1w = (cw * t + aw * (1. - t));
		int p2x = (int)(dx * t + bx * (1. - t) + 0.5);
		double p2w = (dw * t + bw * (1. - t));

		assert(p1x <= p2x);
		draw_horizontal_line(p1x, p1w, p2x, p2w, y, color);
	}
	return true;
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_horizontal_line(int ax, double aw, int bx, double bw, int y, int color)
{
	assert(ax <= bx);

	if (bx < 0)
		return;
	if (ax >= _Xmax)
		return;

	if (y < 0)
		return;
	if (y >= _Ymax)
		return;

	// cut the line if out of screen
	if (ax < 0)
	{
		double t = (double)(bx - 0) / (bx - ax);
		assert(t >= 0.);
		assert(t <= 1.);

		ax = 0; // bx* t + (1. - t) * ax;
		aw = bw * t + (1. - t) * aw;
	}
	if (bx >= _Xmax)
	{
		double t = (float)(bx - (_Xmax - 1)) / (bx - ax);
		assert(t >= 0.);
		assert(t <= 1.);

		bx = _Xmax - 1; // ax* t + (1. - t) * bx;
		bw = aw * t + (1. - t) * bw;
	}

	int* pLine = _pixelBuffer + y * _Xmax;
	float* pWBuffer = _wBuffer + y * _Xmax;

	for (int i = ax; i <= bx; i++)
	{
		float t;
		if (ax == bx)
			t = 1.;
		else
			t = (float)(i - ax) / (bx - ax);

		assert(t >= 0.f);
		assert(t <= 1.f);

		float w = (float)(bw * t + aw * (1.f - t));  //todo optimize

		if (pWBuffer[i] < w)
		{
			pWBuffer[i] = w;
			pLine[i] = color;
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::clear()
{
	for (int j = 0; j < _Ymax; j++)
	{
		int* pixels = j * _Xmax + _pixelBuffer;
		float* w = j * _Xmax + _wBuffer;
		for (int i = 0; i < _Xmax; i++)
		{
			pixels[i] = _iBackgroundColor;
			w[i] = 0.; // z=inf -> w=0
		}
	}
}
void Renderer::set_background(int iBackgroundColor)
{
	_iBackgroundColor = iBackgroundColor;
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_line(const Point3& p1, const Point3& p2, int color)
{
	//todo optimize
	int x1, y1, x2, y2;
	double fx, fy, dx, dy;
	int i, im, decal;
	double zp1, zp2, fz, dz;
	double wEpsilon = 1.e-5; // todo exact formula

	_camera.project(p1, x1, y1, zp1);
	_camera.project(p2, x2, y2, zp2);

	dx = x2 - x1;
	dy = y2 - y1;
	if (fabs(dx) > fabs(dy))
		im = int(fabs(dx));
	else
		im = int(fabs(dy));

	if (im < 1)
		return;

	fx = x1;
	dx = dx / im;
	fy = y1;
	dy = dy / im;
	fz = zp1;
	dz = (zp2 - zp1) / im;

	for (i = 0; i < im; i++)
	{
		if ((fx >= 0) && (fx < _Xmax) && (fy >= 0) && (fy < _Ymax))// && (fz>0) )
		{
			decal = int(fx) + int(fy) * _Xmax;

			if (_wBuffer[decal] < fz + wEpsilon)
			{
				_wBuffer[decal] = (float)fz;
				_pixelBuffer[decal] = color;
			}
		}
		fx += dx;
		fy += dy;
		fz += dz;
	}
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_pixel(const Point3& pPixels, int col)
{
	int ix, iy;
	double zpx;
	float* zp;

	_camera.project(pPixels, ix, iy, zpx);

	if ((ix < 0) || (ix >= _Xmax) || (iy < 0) || (iy >= _Ymax))
		return;

	//teste if visible in wbuffer
	zp = _wBuffer + (int)ix + (int)iy * _Xmax;

	//compare with wbuffer and write if nearer
	if ((*zp) < zpx)
	{
		*zp = (float)zpx;
		_pixelBuffer[(int)ix + (int)iy * _Xmax] = col;
	}
}
////////////////////////////////////////////////////////////////////////////////
void Renderer::draw_polyline(const vector<Point3>& pl, int color)
{
	if (pl.size() < 2)
		return;

	for (int i = 1; i < pl.size(); i++)
		draw_line(pl[i - 1], pl[i], color);
}
////////////////////////////////////////////////////////////////////////////////