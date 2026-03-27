#include <iostream>
using namespace std;

#include "Image.h"
#include "ImageGenerator.h"
#include "ImageIoBmp.h"
#include "ImageUtil.h"
#include "Transform.h"
#include "MeshUtil.h"
#include "OBJFile.h"

////////////////////////////////////////////////////////////////////////////////
int main()
{
	//create a mesh using generated Mandelbrot set 

	cout << "Generating Mandelbrot set" << endl;
	Image imgMandelbrot;
	ImageGenerator::Mandelbrot(imgMandelbrot, 256, 256,128);
	ImageIoBmp::write("test_mesh_mandelbrot.bmp", &imgMandelbrot);
	Image imBW;
	ImageUtil::convert_to_bw(imgMandelbrot, imBW);

	cout << "Converting Mandelbrot to Mesh" << endl;
	Mesh m;
	vector<double> vd;
	imBW.to_double(vd);
	MeshUtil::create_surface_from_z(vd, 256, 256, m);
	m.apply_transform(Scale(1, 1, 0.3));
	
	cout << "Saving Mesh to Obj file" << endl;
	OBJFile::save("sample_mesh_mandelbrot.obj", m);

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
