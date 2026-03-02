#include <iostream>
using namespace std;

#include "Image.h"
#include "ImageGenerator.h"

#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "MeshUtil.h"
#include "Transform.h"
#include "ImageUtil.h"
#include "OBJFile.h"
#include "StepFile.h"

////////////////////////////////////////////////////////////////////////////////
int main()
{
	//create a nurbs using generated Mandelbrot set, save as Nurbs and Mesh

	cout << "Generating Mandelbrot set" << endl;
	Image imgMandelbrot;
	ImageGenerator::Mandelbrot(imgMandelbrot, 64, 64,255);
	Image imBW;
	ImageUtil::convert_to_bw(imgMandelbrot, imBW);

	cout << "Converting Mandelbrot to NurbsSurface" << endl;
	NurbsSurface n;
	vector<double> vd;
	imBW.to_double(vd);
	NurbsUtil::create_from_z(vd, 64, 64, 3, n);
	n.apply_transform(Scale(1., 1., 20. / 255.));

	cout << "Saving NurbsSurface to Obj and step file" << endl;

	Mesh m;
	NurbsUtil::to_mesh(n, m);
	OBJFile::save("sample_nurbs_mandelbrot.obj", m);

	StepWriter sw;
	sw.open("sample_nurbs_mandelbrot.step");
	sw.write(n);

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
