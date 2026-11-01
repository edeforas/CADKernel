#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "Mesh.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

void test_near(double a, double ref, double epsilon=1.e-10,const string& sMessage="")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_square()
{
	cout << endl << "test_nurbssurface_square" << endl;

	NurbsSurface n;
	int degree = 1;
	vector<Point3> points = { Point3(0.,0.,0.),Point3(1.,0.,0.),Point3(0.,1.,0.),Point3(1.,1.,0.) };

	n.set_degree(degree, degree);
	n.set_points(points,2,2);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 p;
			n.evaluate(u, v, p);
			test_near(p.x(), u, 1e-10);
			test_near(p.y(), v, 1e-10);
			test_near(p.z(), 0., 1e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_adaptive_tessellation()
{
	NurbsSurface n;
	n.set_degree(2, 1);
	n.set_points({
		Point3(-1., 0., 0.), Point3(0., 0., 2.), Point3(1., 0., 0.),
		Point3(-1., 1., 0.), Point3(0., 1., 2.), Point3(1., 1., 0.)
	}, 3, 2);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();

	Mesh adaptive;
	NurbsUtil::to_mesh(n, adaptive);
	Mesh coarse;
	NurbsUtil::to_mesh(n, coarse, 1);
	test_bool(adaptive.nb_triangles() > coarse.nb_triangles(),
		"adaptive tessellation should refine a curved surface");
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_differential_plane()
{
	cout << endl << "test_nurbssurface_differential_plane" << endl;

	NurbsSurface n;
	n.set_degree(1, 1);
	n.set_points({ Point3(0.,0.,0.),Point3(1.,0.,0.),Point3(0.,1.,0.),Point3(1.,1.,0.) }, 2, 2);
	n.set_uniform_u();
	n.set_uniform_v();
	n.set_equals_weights();

	Point3 du, dv, duu, duv, dvv;
	n.evaluate_derivatives(0.3, 0.7, du, dv, duu, duv, dvv);
	test_near(du.x(), 1., 1.e-10);
	test_near(du.y(), 0., 1.e-10);
	test_near(dv.x(), 0., 1.e-10);
	test_near(dv.y(), 1., 1.e-10);
	test_near(duu.norm(), 0., 1.e-10);
	test_near(duv.norm(), 0., 1.e-10);
	test_near(dvv.norm(), 0., 1.e-10);

	Point3 normal;
	test_bool(n.normal(0.4, 0.6, normal));
	test_near(std::fabs(normal.z()), 1., 1.e-10);

	double K, H;
	test_bool(n.curvature(0.4, 0.6, K, H));
	test_near(K, 0., 1.e-10);
	test_near(H, 0., 1.e-10);
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_differential_cylinder()
{
	cout << endl << "test_nurbssurface_differential_cylinder" << endl;

	NurbsSurface n;
	n.set_degree(2, 1);
	n.set_knots_u({ 0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1. });
	n.set_knots_v({ 0., 0., 1., 1. });
	n.set_weights({
		1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,
		1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.
	});
	n.set_points({
		Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.),
		Point3(-1.,1.,0.),Point3(-1.,0.,0.),Point3(-1.,-1.,0.),
		Point3(0.,-1.,0.),Point3(1.,-1.,0.),Point3(1.,0.,0.),
		Point3(1.,0.,2.),Point3(1.,1.,2.),Point3(0.,1.,2.),
		Point3(-1.,1.,2.),Point3(-1.,0.,2.),Point3(-1.,-1.,2.),
		Point3(0.,-1.,2.),Point3(1.,-1.,2.),Point3(1.,0.,2.)
	}, 9, 2);

	for (double u = 0.1; u <= 0.9; u += 0.2)
	{
		Point3 normal;
		test_bool(n.normal(u, 0.5, normal));
		test_near(normal.z(), 0., 1.e-8);

		double K, H;
		test_bool(n.curvature(u, 0.5, K, H));
		test_near(K, 0., 2.e-5);
		test_near(std::fabs(H), 0.5, 2.e-3);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_elevate_degree_from_deg1()
{
	// test all degree elevation combinations up to 3
	cout << endl << "test_nurbssurface_elevate_degree_from_deg1" << endl;

	int nbPointsU = 11;
	int nbPointsV = 11;

	NurbsSurface ndeg_u1v1;
	int degree = 1;
	vector<Point3> points;
	vector<double> weights;

	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
		{
			points.push_back(Point3(u, v, (float)rand() / RAND_MAX * 3.));
			weights.push_back((float)rand() / RAND_MAX + 1);
		}

	ndeg_u1v1.set_degree(degree, degree);
	ndeg_u1v1.set_points(points, nbPointsU, nbPointsV);
	ndeg_u1v1.set_weights(weights);

	vector<double> knot_u = {0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1};
	vector<double> knot_v = { 0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1};
	ndeg_u1v1.set_knots_u(knot_u);
	ndeg_u1v1.set_knots_v(knot_v);

	NurbsSurface ndeg_u2v1 = ndeg_u1v1;
	ndeg_u2v1.degree_elevation_u();

	NurbsSurface ndeg_u3v1 = ndeg_u2v1;
	ndeg_u3v1.degree_elevation_u();

	NurbsSurface ndeg_u1v2 = ndeg_u1v1;
	ndeg_u1v2.degree_elevation_v();

	NurbsSurface ndeg_u1v3 = ndeg_u1v2;
	ndeg_u1v3.degree_elevation_v();

	NurbsSurface ndeg_u2v2 = ndeg_u2v1;
	ndeg_u2v2.degree_elevation_v();

	NurbsSurface ndeg_u3v2 = ndeg_u2v2;
	ndeg_u3v2.degree_elevation_u();

	NurbsSurface ndeg_u2v3 = ndeg_u2v2;
	ndeg_u2v3.degree_elevation_v();

	NurbsSurface ndeg_u3v3 = ndeg_u2v3;
	ndeg_u3v3.degree_elevation_u();

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 pdeg_u1v1,pdeg_u2v1, pdeg_u3v1, pdeg_u1v2, pdeg_u1v3, pdeg_u2v2, pdeg_u2v3, pdeg_u3v2, pdeg_u3v3;
			ndeg_u1v1.evaluate(u, v, pdeg_u1v1);

			ndeg_u2v1.evaluate(u, v, pdeg_u2v1);
			ndeg_u3v1.evaluate(u, v, pdeg_u3v1);
			ndeg_u1v2.evaluate(u, v, pdeg_u1v2);
			ndeg_u1v3.evaluate(u, v, pdeg_u1v3);
			ndeg_u2v2.evaluate(u, v, pdeg_u2v2);
			ndeg_u2v3.evaluate(u, v, pdeg_u2v3);
			ndeg_u3v2.evaluate(u, v, pdeg_u3v2);
			ndeg_u3v3.evaluate(u, v, pdeg_u3v3);

			test_near((pdeg_u1v1 - pdeg_u2v1).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u3v1).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u1v2).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u1v3).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u2v2).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u2v3).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u3v2).norm(), 0, 1.e-10);
			test_near((pdeg_u1v1 - pdeg_u3v3).norm(), 0, 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_elevate_degree_from_deg2()
{
	// test all degree elevation combinations up to 3
	cout << endl << "test_nurbssurface_elevate_degree_from_deg2" << endl;

	int nbPointsU = 11;
	int nbPointsV = 11;

	NurbsSurface ndeg_u2v2;
	int degree = 2;
	vector<Point3> points;
	vector<double> weights;

	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
		{
			points.push_back(Point3(u, v, (float)rand() / RAND_MAX * 3.));
			weights.push_back((float)rand() / RAND_MAX + 1);
		}

	ndeg_u2v2.set_degree(degree, degree);
	ndeg_u2v2.set_points(points, nbPointsU, nbPointsV);
	ndeg_u2v2.set_weights(weights);

	vector<double> knot_u = { 0,0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1 };
	vector<double> knot_v = { 0,0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1 };
	ndeg_u2v2.set_knots_u(knot_u);
	ndeg_u2v2.set_knots_v(knot_v);

	NurbsSurface ndeg_u3v2 = ndeg_u2v2;
	ndeg_u3v2.degree_elevation_u();

	NurbsSurface ndeg_u2v3 = ndeg_u2v2;
	ndeg_u2v3.degree_elevation_v();

	NurbsSurface ndeg_u3v3 = ndeg_u2v3;
	ndeg_u3v3.degree_elevation_u();

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 pdeg_u2v2, pdeg_u2v3, pdeg_u3v2, pdeg_u3v3;

			ndeg_u2v2.evaluate(u, v, pdeg_u2v2);
			ndeg_u2v3.evaluate(u, v, pdeg_u2v3);
			ndeg_u3v2.evaluate(u, v, pdeg_u3v2);
			ndeg_u3v3.evaluate(u, v, pdeg_u3v3);

			test_near((pdeg_u2v2 - pdeg_u2v3).norm(), 0, 1.e-10);
			test_near((pdeg_u2v2 - pdeg_u3v2).norm(), 0, 1.e-10);
			test_near((pdeg_u2v2 - pdeg_u3v3).norm(), 0, 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_deg1_insert_knot()
{
	cout << endl << "test_nurbssurface_deg1_insert_knot" << endl;

	int nbPointsU = 11;
	int nbPointsV = 11;

	NurbsSurface ndeg_u1v1;
	int degree = 1;
	vector<Point3> points;
	vector<double> weights;

	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
		{
			points.push_back(Point3(u, v, (float)rand() / RAND_MAX * 3.));
			weights.push_back((float)rand() / RAND_MAX + 1);
		}

	ndeg_u1v1.set_degree(degree, degree);
	ndeg_u1v1.set_points(points, nbPointsU, nbPointsV);
	ndeg_u1v1.set_weights(weights);

	vector<double> knot_u = { 0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1 };
	vector<double> knot_v = { 0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1 };
	ndeg_u1v1.set_knots_u(knot_u);
	ndeg_u1v1.set_knots_v(knot_v);

	NurbsSurface ndeg_u1v1_knot = ndeg_u1v1;
	ndeg_u1v1_knot.insert_knot_u(0.15);
	ndeg_u1v1_knot.insert_knot_u(0.45);
	ndeg_u1v1_knot.insert_knot_u(0.88);
	ndeg_u1v1_knot.insert_knot_v(0.25);
	ndeg_u1v1_knot.insert_knot_v(0.66);
	ndeg_u1v1_knot.insert_knot_v(0.77);

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 pdeg_u1v1, pdeg_u1v1_knot;
			ndeg_u1v1.evaluate(u, v, pdeg_u1v1);
			ndeg_u1v1_knot.evaluate(u, v, pdeg_u1v1_knot);

			test_near((pdeg_u1v1 - pdeg_u1v1_knot).norm(), 0, 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_deg2_insert_knot()
{
	cout << endl << "test_nurbssurface_deg2_insert_knot" << endl;

	int nbPointsU = 11;
	int nbPointsV = 11;

	NurbsSurface ndeg_u2v2;
	int degree = 2;
	vector<Point3> points;
	vector<double> weights;

	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
		{
			points.push_back(Point3(u, v, (float)rand() / RAND_MAX * 3.));
			weights.push_back((float)rand() / RAND_MAX + 1);
		}

	ndeg_u2v2.set_degree(degree, degree);
	ndeg_u2v2.set_points(points, nbPointsU, nbPointsV);
	ndeg_u2v2.set_weights(weights);

	vector<double> knot_u = { 0,0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1 };
	vector<double> knot_v = { 0,0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1 };
	ndeg_u2v2.set_knots_u(knot_u);
	ndeg_u2v2.set_knots_v(knot_v);

	NurbsSurface ndeg_u2v2_knot = ndeg_u2v2;
	ndeg_u2v2_knot.insert_knot_u(0.15);
	ndeg_u2v2_knot.insert_knot_u(0.45);
	ndeg_u2v2_knot.insert_knot_u(0.88);
	ndeg_u2v2_knot.insert_knot_v(0.25);
	ndeg_u2v2_knot.insert_knot_v(0.66);
	ndeg_u2v2_knot.insert_knot_v(0.77);

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 pdeg_u1v1, pdeg_u1v1_knot;
			ndeg_u2v2.evaluate(u, v, pdeg_u1v1);
			ndeg_u2v2_knot.evaluate(u, v, pdeg_u1v1_knot);

			test_near((pdeg_u1v1 - pdeg_u1v1_knot).norm(), 0, 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_deg3_insert_knot()
{
	cout << endl << "test_nurbssurface_deg3_insert_knot" << endl;

	int nbPointsU = 11;
	int nbPointsV = 11;

	NurbsSurface ndeg_u3v3;
	int degree = 3;
	vector<Point3> points;
	vector<double> weights;

	for (int v = 0; v < nbPointsV; v++)
		for (int u = 0; u < nbPointsU; u++)
		{
			points.push_back(Point3(u, v, (float)rand() / RAND_MAX * 3.));
			weights.push_back((float)rand() / RAND_MAX + 1);
		}

	ndeg_u3v3.set_degree(degree, degree);
	ndeg_u3v3.set_points(points, nbPointsU, nbPointsV);
	ndeg_u3v3.set_weights(weights);

	vector<double> knot_u = { 0,0,0,0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1,1 };
	vector<double> knot_v = { 0,0,0,0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1,1 };
	ndeg_u3v3.set_knots_u(knot_u);
	ndeg_u3v3.set_knots_v(knot_v);

	NurbsSurface ndeg_u3v3_knot = ndeg_u3v3;
	ndeg_u3v3_knot.insert_knot_u(0.15);
	ndeg_u3v3_knot.insert_knot_u(0.45);
	ndeg_u3v3_knot.insert_knot_u(0.88);
	ndeg_u3v3_knot.insert_knot_v(0.25);
	ndeg_u3v3_knot.insert_knot_v(0.66);
	ndeg_u3v3_knot.insert_knot_v(0.77);

	for (double u = 0.; u <= 1.; u += 0.1)
		for (double v = 0.; v <= 1.; v += 0.1)
		{
			Point3 pdeg_u1v1, pdeg_u1v1_knot;
			ndeg_u3v3.evaluate(u, v, pdeg_u1v1);
			ndeg_u3v3_knot.evaluate(u, v, pdeg_u1v1_knot);

			test_near((pdeg_u1v1 - pdeg_u1v1_knot).norm(), 0, 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
void test_nurbssurface_assignment_mixed_degrees()
{
	cout << endl << "test_nurbssurface_assignment_mixed_degrees" << endl;

	NurbsSurface src;
	src.set_degree(2, 1);
	src.set_knots_u({ 0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1. });
	src.set_knots_v({ 0., 0., 1., 1. });
	src.set_weights({
		1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,
		1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.
	});
	src.set_points({
		Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.),
		Point3(-1.,1.,0.),Point3(-1.,0.,0.),Point3(-1.,-1.,0.),
		Point3(0.,-1.,0.),Point3(1.,-1.,0.),Point3(1.,0.,0.),
		Point3(1.,0.,2.),Point3(1.,1.,2.),Point3(0.,1.,2.),
		Point3(-1.,1.,2.),Point3(-1.,0.,2.),Point3(-1.,-1.,2.),
		Point3(0.,-1.,2.),Point3(1.,-1.,2.),Point3(1.,0.,2.)
	}, 9, 2);

	NurbsSurface dst;
	dst.set_degree(1, 3);
	dst = src;

	for (double u = 0.05; u <= 0.95; u += 0.15)
		for (double v = 0.05; v <= 0.95; v += 0.15)
		{
			Point3 pSrc, pDst;
			src.evaluate(u, v, pSrc);
			dst.evaluate(u, v, pDst);
			test_near((pSrc - pDst).norm(), 0., 1.e-10);
		}
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	test_nurbssurface_square();
	test_nurbssurface_adaptive_tessellation();
	test_nurbssurface_differential_plane();
	test_nurbssurface_differential_cylinder();

	test_nurbssurface_elevate_degree_from_deg1();
	test_nurbssurface_elevate_degree_from_deg2();
	
	test_nurbssurface_deg1_insert_knot();
	test_nurbssurface_deg2_insert_knot();
	test_nurbssurface_deg3_insert_knot();
	test_nurbssurface_assignment_mixed_degrees();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////