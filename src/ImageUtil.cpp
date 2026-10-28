#include "ImageUtil.h"

#include "Image.h"
#include "Mesh.h"

#include <cassert>
#include <algorithm>
#include <cmath>
using namespace std;

///////////////////////////////////////////////////////////////////////////
void ImageUtil::convert_to_grayscale(const Image& im, Image& imBW)
{
	// todo check im format before conversion

	imBW.resize(im.width(), im.height(), 1);

	for (int i=0; i< im.width();i++)
		for (int j = 0; j < im.height(); j++)
		{
			imBW(i, j, 0) = (im(i, j, 0) + im(i, j, 1) + im(i, j, 2)) / 3;
		}
}
///////////////////////////////////////////////////////////////////////////
void ImageUtil::scale(Image& im, double dScale)
{
	for (int i = 0; i < im.width(); i++)
		for (int j = 0; j < im.height(); j++)
			for (int p= 0;p< im.planes();p++)
				im(i, j, p) = (unsigned char)(im(i, j, p) * dScale + 0.5);
}
///////////////////////////////////////////////////////////////////////////
