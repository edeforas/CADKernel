#ifndef _Renderer_
#define _Renderer_

#include "RendererCamera.h"
#include "RendererLight.h"
#include "Geometry.h"
class Mesh;
class NurbsSurface;
class NurbsSolid;

#include <vector>
using namespace std;

class Renderer
{
public:
	Renderer(int* pBuffer, int width, int height);
	~Renderer();

	void clear();
	void set_background(int iBackgroundColor);
	void set_camera(double ox, double oy, double oz, double ahead, double yaw, double pitch, double roll, double zoom);
	void add_ambient_light(int iAmbiantColor, double dAmbiantFactor = 1.);
	void add_diffuse_light(int iDiffuseColor, double dDiffuseFactor, const Point3& direction);

	void draw_mesh(const Mesh& m, bool bDrawEdges = false, int iColor = -1); //use mesh color if iColor is -1 
	void draw_surface(const NurbsSurface& n, int iNbSegments = 5, bool bDrawEdges = false, int iColor = -1);
	void draw_solid(const NurbsSolid& n, int iNbSegments = 5, bool bDrawEdges = false, int iColor = -1);
	bool draw_triangle_1color(const Point3& A, const Point3& B, const Point3& C, int color, bool bTwofaces = false); //return true if face was visible
	void draw_line(const Point3& p1, const Point3& p2, int color);
	void draw_pixel(const Point3& pPixels, int color);
	void draw_polyline(const vector<Point3>& pl, int color);

private:
	bool draw_trapeze(double ax, double aw, double bx, double bw, int ay, double cx, double cw, double dx, double dw, int cy, int color); //return true if face was visible
	void draw_horizontal_line(int ax, double aw, int bx, double bw, int y, int color);

	int compute_color_with_lights(int iColor, const Point3& normal);

	RendererCamera _camera;
	vector<RendererLight*> _lights;

	int* _pixelBuffer;
	float* _wBuffer;
	int _Xmax, _Ymax;
	int _iBackgroundColor;
};

#endif
