#include <iostream>
using namespace std;

#include "Image.h"
#include "ImageGenerator.h"

#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "NurbsFactory.h"
#include "NurbsRuled.h"
#include "NurbsSolid.h"
#include "Transform.h"
#include "ImageUtil.h"
#include "OBJFile.h"
#include "StepFile.h"

////////////////////////////////////////////////////////////////////////////////
int main()
{
	//create a Solid 3D Mandelbrot set, save as Nurbs and Mesh

	cout << "Generate Mandelbrot set" << endl;
	Image imgMandelbrot;
	ImageGenerator::Mandelbrot(imgMandelbrot, 64, 64,255);
	Image imBW;
	ImageUtil::convert_to_bw(imgMandelbrot, imBW);

	cout << "Convert Mandelbrot to NurbsSurface" << endl;
	NurbsSurface n;
	vector<double> vd;
	imBW.to_double(vd);
	NurbsUtil::create_surface_from_z(vd, 64, 64, 3, n);
	n.apply_transform(Scale(1., 1., 20. / 255.));

	cout << "Create base" << endl;
	NurbsSurface nsBase;
	NurbsFactory::create_quad(Point3(0, 0, -10), Point3(0, 63, -10), Point3(63, 63, -10), Point3(63, 0, -10), nsBase);

	cout << "Create edges" << endl;
	NurbsRuled nr;
	NurbsSurface nsU0, nsU1, nsV0, nsV1;
	nr.create_ruled(n.edge_u0(), nsBase.edge_u0(),nsU0);
	nr.create_ruled(n.edge_u1(), nsBase.edge_u1(), nsU1);
	nr.create_ruled(n.edge_v0(), nsBase.edge_v0(), nsV0);
	nr.create_ruled(n.edge_v1(), nsBase.edge_v1(), nsV1);

	cout << "Assemble surfaces in a solid" << endl;
	NurbsSolid ns;
	ns.add_surface(n);
	ns.add_surface(nsBase);
	ns.add_surface(nsU0);
	ns.add_surface(nsU1);
	ns.add_surface(nsV0);
	ns.add_surface(nsV1);

	cout << "Save in .step and .obj format" << endl;

	Mesh m;
	NurbsUtil::to_mesh(ns, m);
	OBJFile::save("sample_nurbs_mandelbrot.obj", m);

	StepWriter sw;
	sw.open("sample_nurbs_mandelbrot.step");
	sw.write(ns);
	sw.close();

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
