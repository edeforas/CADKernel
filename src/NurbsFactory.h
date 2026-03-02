#ifndef NurbsFactory_
#define NurbsFactory_

#include "Geometry.h"

class NurbsCurve; 
class NurbsSurface;
class NurbsSolid;

#include <vector>

///////////////////////////////////////////////////////////////////////////
namespace NurbsFactory
{
	// curves
	void create_circle(double dRadius, NurbsCurve& nc);
	void create_filled_surface(const NurbsCurve& nc, NurbsSurface& ns);
	
	// surfaces
	void create_triangle(const Point3& p1, const Point3& p2, const Point3& p3, NurbsSurface& ns);
	void create_quad(const Point3& p1, const Point3& p2, const Point3& p3, const Point3& p4, NurbsSurface& ns);
	void create_disk(double dRadius, NurbsSurface& ns);
	void create_sphere(double dRadius, NurbsSurface& ns);
	void create_torus(double dMajorRadius, double dMinorRadius, NurbsSurface& ns);

	//solid
	void create_box(double dSizeX, double dSizeY, double dSizeZ, NurbsSolid& n);
	void create_cylinder(double dRadius,double dHeight, NurbsSolid& n);
	void create_sphere(double dRadius, NurbsSolid& n);
	void create_torus(double dMajorRadius, double dMinorRadius, NurbsSolid& n);
}
///////////////////////////////////////////////////////////////////////////

#endif