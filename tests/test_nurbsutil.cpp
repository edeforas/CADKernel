#include "NurbsCurve.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "StepFile.h"

#include "MeshFactory.h"
#include "OBJFile.h"

#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;

void test_near(double a, double ref, double epsilon=1.e-10,const string& sMessage="")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsutil_create_icosahedron()
{
	Mesh ico,ico_out;
	MeshFactory::create_icosahedron(30.,ico);
	
	NurbsSolid n;
	NurbsUtil::create_from_mesh(ico,n);

	NurbsUtil::to_mesh(n, ico_out, 10);
	OBJFile::save("test_nurbsutil_create_icosahedron_orig.obj", ico);
	OBJFile::save("test_nurbsutil_create_icosahedron_out.obj", ico_out);

	StepWriter sw;
	sw.open("test_nurbsutil_create_icosahedron.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsutil_create_dodecahedron()
{
	Mesh dodec, dodec_out;
	MeshFactory::create_dodecahedron(30., dodec);

	NurbsSolid n;
	NurbsUtil::create_from_mesh(dodec, n);

	NurbsUtil::to_mesh(n, dodec_out, 10);
	OBJFile::save("test_nurbsutil_create_dodecahedron_orig.obj", dodec);
	OBJFile::save("test_nurbsutil_create_dodecahedron_out.obj", dodec_out);

	StepWriter sw;
	sw.open("test_nurbsutil_create_dodecahedron.step");
	sw.write(n);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbsutil_create_curve_from_points_any_degree()
{
	cout << endl << "test_nurbsutil_create_curve_from_points_any_degree" << endl;

	vector<Point3> points = {
		Point3(0., 0., 0.),
		Point3(2., 1., 0.),
		Point3(4., 0., 0.),
		Point3(6., 2., 0.),
		Point3(8., 0., 0.)
	};

	NurbsCurve n;
	NurbsUtil::create_curve_from_points(points, 3, n);
	test_near((double)n.degree(), 3., 1.e-12, "degree should be set to requested value when valid");
	test_near((double)n.nb_points(), (double)points.size(), 1.e-12, "curve should keep all points");

	Point3 pStart, pEnd;
	n.evaluate(0., pStart);
	n.evaluate(1., pEnd);
	test_near((pStart - points.front()).norm(), 0., 1.e-8, "curve start should match first control point");
	test_near((pEnd - points.back()).norm(), 0., 1.e-8, "curve end should match last control point");

	NurbsUtil::create_curve_from_points(points, 100, n);
	test_near((double)n.degree(), (double)points.size() - 1., 1.e-12, "degree should be clamped to points.size()-1");
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	test_nurbsutil_create_icosahedron();
	test_nurbsutil_create_dodecahedron();
	test_nurbsutil_create_curve_from_points_any_degree();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////