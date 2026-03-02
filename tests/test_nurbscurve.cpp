#include "NurbsCurve.h"
#include "NurbsUtil.h"
#include "NurbsFactory.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
using namespace std;

void test(bool b, string s = "")
{
	if (b == false)
	{
		cout << "condition: [" << s << "] is not realized!" << endl;
		exit(-1);
	}
}

void test_near(double a, double ref, double epsilon = 1.e-10, const string& sMessage = "")
{
	if ((a > ref + epsilon) || (a < ref - epsilon))
	{
		cerr << "Test Error: " << sMessage.c_str() << "value=" << a << " ref=" << ref << endl;
		exit(-1);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_knots_span()
{
	cout << endl << "test_knots_span" << endl;

	NurbsCurve n;
	vector<double> knots = { 0,0,1,1 };

	int i_0 = n.find_knot_span(knots, 0.);
	test_near(i_0, 1);

	int i_05 = n.find_knot_span(knots, 0.5);
	test_near(i_05, 1);

	int i_1 = n.find_knot_span(knots, 1.);
	test_near(i_1, 1);

	//extrapolation
	int i_m01 = n.find_knot_span(knots, -0.1);
	test_near(i_m01, 1);

	int i_11 = n.find_knot_span(knots, 1.1);
	test_near(i_11, 1);
}
///////////////////////////////////////////////////////////////////////////
void test_simple_segment()
{
	cout << endl << "test_simple_segment" << endl;

	NurbsCurve n;
	int degree = 1;
	vector<double> knots = { 0,0,1,1 };
	vector<double> weights(2, 1.);
	vector<Point3> points = { Point3(0.,0.,0.),Point3(1.,1.,0.) };

	n.set_degree(degree);
	n.set_knots(knots);
	n.set_weights(weights);
	n.set_points(points);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		test_near(p.x(), u, 1e-10);
		test_near(p.y(), u, 1e-10);
		test_near(p.z(), 0, 1e-10);

		//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << endl;
	}
}
///////////////////////////////////////////////////////////////////////////
void test_curve_derivatives_line()
{
	cout << endl << "test_curve_derivatives_line" << endl;

	NurbsCurve n;
	n.set_degree(1);
	n.set_knots({ 0.,0.,1.,1. });
	n.set_weights({ 1.,1. });
	n.set_points({ Point3(0.,0.,0.), Point3(1.,1.,0.) });

	Point3 d1, d2;
	n.evaluate_derivatives(0.3, d1, d2);
	test_near(d1.x(), 1., 1.e-10);
	test_near(d1.y(), 1., 1.e-10);
	test_near(d1.z(), 0., 1.e-10);
	test_near(d2.norm(), 0., 1.e-10);
	test_near(n.curvature(0.5), 0., 1.e-10);

	Point3 tangent, normal;
	test(n.tangent(0.5, tangent));
	test(!n.normal(0.5, normal));
}
///////////////////////////////////////////////////////////////////////////
void test_quarter_circle_deg2()
{
	cout << endl << "test_quarter_circle_deg2" << endl;

	NurbsCurve n;
	int degree = 2;
	vector<double> knots = { 0,0,0,1,1,1 };
	vector<double> weights = { 1.,1. / sqrt(2.),1. };
	vector<Point3> points = { Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.) };

	n.set_degree(degree);
	n.set_knots(knots);
	n.set_weights(weights);
	n.set_points(points);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		double norm = p.norm();
		//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
		test_near(norm, 1., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_curve_differential_circle()
{
	cout << endl << "test_curve_differential_circle" << endl;

	NurbsCurve n;
	n.set_degree(2);
	n.set_knots({ 0.,0.,0.,1.,1.,1. });
	n.set_weights({ 1.,1. / sqrt(2.),1. });
	n.set_points({ Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.) });

	for (double u = 0.1; u <= 0.9; u += 0.2)
	{
		Point3 p;
		n.evaluate(u, p);

		Point3 t, nn;
		test(n.tangent(u, t), "tangent should be defined on circle arc");
		test(n.normal(u, nn), "normal should be defined on circle arc");

		test_near(t.norm(), 1., 1.e-8);
		test_near(nn.norm(), 1., 1.e-8);
		test_near(t.dot_product(p), 0., 1.e-8);
		test_near(n.curvature(u), 1., 1.e-6);

		const Point3 r = p.normalized();
		test_near(nn.dot_product(r), -1., 1.e-5);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_half_circle_deg3()
{
	cout << endl << "test_half_circle_deg3" << endl;

	NurbsCurve n;
	int degree = 3;
	vector<double> knots = { 0,0,0,0,1,1,1,1 };
	vector<double> weights = { 1.,1. / 3.,1. / 3.,1. };
	vector<Point3> points = { Point3(1.,0.,0.),Point3(1.,2.,0.),Point3(-1.,2.,0.),Point3(-1.,0.,0.) };

	n.set_degree(degree);
	n.set_knots(knots);
	n.set_weights(weights);
	n.set_points(points);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		double norm = p.norm();
		//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
		test_near(norm, 1., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_circle_approximation_deg3()
{
	//	approximation using cubic NURBS, not rational
	//	from "Good Approximation of Circles by Curvature-Continuous Bezier curves"
	//	0.55 number optimized to 0.5519151 for a lower max error

	cout << endl << "test_circle_approximation_deg3" << endl;

	NurbsCurve n;
	int degree = 3;
	vector<double> knots = { 0,0,0,0,1,1,1,1 };
	vector<double> weights(4, 1.);
	vector<Point3> points = { Point3(1.,0.,0.),Point3(1.,0.5519151,0.), Point3(0.5519151,1.,0.), Point3(0.,1.,0.) };

	n.set_degree(degree);
	n.set_knots(knots);
	n.set_weights(weights);
	n.set_points(points);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		double norm = p.norm();
		//cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
		test_near(norm, 1., 2.e-4);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_full_circle_deg2()
{
	cout << endl << "test_full_circle_deg2" << endl;

	NurbsCurve n;
	int degree = 2;
	vector<double> knots = { 0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1. };
	vector<double> weights = { 1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1.,1. / sqrt(2.),1. };
	vector<Point3> points = {
		Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.),
		Point3(-1.,1.,0.),Point3(-1.,0.,0.),Point3(-1.,-1.,0.),
		Point3(0.,-1.,0.),Point3(1.,-1.,0.),Point3(1.,0.,0.)
	};

	n.set_degree(degree);
	n.set_knots(knots);
	n.set_weights(weights);
	n.set_points(points);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		double norm = p.norm();
		//		cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << " norm=" << norm << endl;
		test_near(norm, 1., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_random_deg2()
{
	cout << endl << "test_random_deg2" << endl;

	NurbsCurve n;
	int degree = 2;
	vector<Point3> points;
	int nbPoints = 20;

	for (int i = 0; i < nbPoints; i++)
		points.push_back(Point3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX));

	NurbsUtil::create_curve_from_points(points, degree, n);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p;
		n.evaluate(u, p);
		//	cout << "u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << endl;
	}
}
///////////////////////////////////////////////////////////////////////////
void test_knot_insertion_deg1()
{
	cout << endl << "test_knot_insertion_deg1" << endl;

	NurbsCurve n_deg1;
	int degree = 1;
	vector<Point3> points = { Point3(2.,7.,0.1),Point3(1.,-4.,8.),Point3(-2.,0.,0.6) };
	vector<double> weights = { 1.,1.2,0.7 };
	vector<double> knots = { 0,0, 0.6 ,1,1 };
	n_deg1.set_degree(degree);
	n_deg1.set_points(points);
	n_deg1.set_weights(weights);
	n_deg1.set_knots(knots);

	// insert knots and compare
	NurbsCurve n_deg1k = n_deg1;
	n_deg1k.insert_knot(0.2);
	n_deg1k.insert_knot(0.7);
	n_deg1k.insert_knot(0.9);

	for (double u = 0.; u <= 1.; u += 0.1)
	{
		Point3 p1, p2;
		n_deg1.evaluate(u, p1);
		n_deg1k.evaluate(u, p2);

		test_near((p1 - p2).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_knot_insertion_deg2()
{
	cout << endl << "test_knot_insertion_deg2" << endl;

	NurbsCurve n_deg2;
	int degree = 2;
	vector<Point3> points = { Point3(2.,7.,0.1),Point3(1.,-4.,8.),Point3(3.,-2.,7.),Point3(-2.,0.,0.6) };
	vector<double> weights = { 1.,1.2,0.9,0.7 };
	vector<double> knots = { 0,0,0, 0.6, 1,1,1 };
	n_deg2.set_degree(degree);
	n_deg2.set_points(points);
	n_deg2.set_weights(weights);
	n_deg2.set_knots(knots);

	// insert knots and compare
	NurbsCurve n_deg2k = n_deg2;
	n_deg2k.insert_knot(0.2);
	n_deg2k.insert_knot(0.7);
	n_deg2k.insert_knot(0.8);

	for (double u = 0.; u <= 1.; u += 0.1)
	{
		Point3 p1, p2;
		n_deg2.evaluate(u, p1);
		n_deg2k.evaluate(u, p2);

		test_near((p1 - p2).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_knot_insertion_deg3()
{
	cout << endl << "test_knot_insertion_deg3" << endl;

	NurbsCurve n_deg3;
	int degree = 3;
	vector<Point3> points = { Point3(2.,7.,0.1),Point3(1.,-4.,8.),Point3(3.,-2.,7.),Point3(4.,-6.,-5.),Point3(-2.,0.,0.6) };
	vector<double> weights = { 1.,1.2,0.9,0.5,0.7 };
	vector<double> knots = { 0,0,0,0, 0.6 , 1,1,1,1 };
	n_deg3.set_degree(degree);
	n_deg3.set_points(points);
	n_deg3.set_weights(weights);
	n_deg3.set_knots(knots);

	// insert knots and compare
	NurbsCurve n_deg3k = n_deg3;
	n_deg3k.insert_knot(0.2);
	n_deg3k.insert_knot(0.7);
	n_deg3k.insert_knot(0.8);

	for (double u = 0.; u <= 1.; u += 0.1)
	{
		Point3 p1, p2;
		n_deg3.evaluate(u, p1);
		n_deg3k.evaluate(u, p2);

		test_near((p1 - p2).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_degree_elevation_from_deg1()
{
	cout << endl << "test_degree_elevation_from_deg1" << endl;

	NurbsCurve n_deg1;
	int degree = 1;

	vector<double> knots = { 0,0,0.6,0.8,1,1 };
	vector<double> weights = { 1.,0.5,0.7,2 };
	vector<Point3> points = {
		Point3(-2,-3,1.),
		Point3(5,2,3.),
		Point3(-5,3,-1.),
		Point3(-1,-4,2.),
	};

	n_deg1.set_degree(degree);
	n_deg1.set_knots(knots);
	n_deg1.set_points(points);
	n_deg1.set_weights(weights);

	// elevate degree from 1 to 2
	NurbsCurve n_deg2 = n_deg1;
	n_deg2.degree_elevation();
	test(n_deg2.degree() == 2);

	// elevate degree from 2 to 3
	NurbsCurve n_deg3 = n_deg2;
	n_deg3.degree_elevation();
	test(n_deg3.degree() == 3);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p1, p2, p3;
		n_deg1.evaluate(u, p1);
		n_deg2.evaluate(u, p2);
		n_deg3.evaluate(u, p3);
		//cout << "p1 u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << endl;
		//cout << "p2 u=" << u << " x=" << p2.x() << " y=" << p2.y() << " z=" << p2.z() << endl;

		test_near((p1 - p2).norm(), 0., 1.e-10);
		test_near((p1 - p3).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_degree_elevation_from_deg2()
{
	cout << endl << "test_degree_elevation_from_deg2" << endl;

	NurbsCurve n_deg2;
	int degree = 2;

	vector<double> knots = { 0,0,0, 0.5, 1,1,1 };
	vector<double> weights = { 1., 0.7 ,0.8, 1.1 };
	vector<Point3> points = {
		Point3(-2,-3,1.),
		Point3(5,2,3.),
		Point3(-5,3,-6.),
		Point3(-1,-4,2.),
	};

	n_deg2.set_degree(degree);
	n_deg2.set_knots(knots);
	n_deg2.set_points(points);
	n_deg2.set_weights(weights);

	// elevate degree from 2 to 3
	NurbsCurve n_deg3 = n_deg2;
	n_deg3.degree_elevation();
	test(n_deg3.degree() == 3);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p2, p3;
		n_deg2.evaluate(u, p2);
		n_deg3.evaluate(u, p3);

		//cout << "p1 u=" << u << " x=" << p.x() << " y=" << p.y() << " z=" << p.z() << endl;
		//cout << "p2 u=" << u << " x=" << p2.x() << " y=" << p2.y() << " z=" << p2.z() << endl;

		test_near((p2 - p3).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_degree_elevation_circle_deg2()
{
	cout << endl << "test_degree_elevation_circle_deg2" << endl;

	NurbsCurve n_deg2;
	NurbsFactory::create_circle(1, n_deg2);
	test(n_deg2.degree() == 2);

	// elevate degree from 2 to 3
	NurbsCurve n_deg3 = n_deg2;
	n_deg3.degree_elevation();
	test(n_deg3.degree() == 3);

	for (double u = 0.; u <= 1.; u += 0.01)
	{
		Point3 p2, p3;
		n_deg2.evaluate(u, p2);
		n_deg3.evaluate(u, p3);

		test_near(p2.norm(), 1., 1.e-10);
		test_near(p3.norm(), 1., 1.e-10);

		test_near((p2 - p3).norm(), 0., 1.e-10);
	}
}
///////////////////////////////////////////////////////////////////////////
void test_is_closed()
{
	cout << endl << "test_is_closed" << endl;

	{
		NurbsCurve n;
		test(n.is_closed() == false);
	}

	{
		NurbsCurve n;
		vector<Point3> points = { Point3(1.,0.,0.) };
		NurbsUtil::create_curve_from_points(points, 1, n);
		test(n.is_closed() == false);
	}

	{
		NurbsCurve n;
		vector<Point3> points = { Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(0.,1.,0.) };
		NurbsUtil::create_curve_from_points(points, 1, n);
		test(n.is_closed() == false);
	}

	{
		NurbsCurve n;
		vector<Point3> points = { Point3(1.,0.,0.),Point3(1.,1.,0.),Point3(1.,0.,0.) };
		NurbsUtil::create_curve_from_points(points, 1, n);
		test(n.is_closed() == true);
	}

	{
		NurbsCurve n;
		vector<Point3> points = { Point3(1. + 1.1e-6,0.,0.),Point3(1.,1.,0.),Point3(1.,0.,0.) };
		NurbsUtil::create_curve_from_points(points, 1, n);
		test(n.is_closed() == false);
	}
}
///////////////////////////////////////////////////////////////////////////
int main()
{
	test_knots_span();
	test_simple_segment();
	test_curve_derivatives_line();
	test_quarter_circle_deg2();
	test_curve_differential_circle();
	test_half_circle_deg3();
	test_full_circle_deg2();
	test_circle_approximation_deg3();

	test_random_deg2();

	test_knot_insertion_deg1();
	test_knot_insertion_deg2();
	test_knot_insertion_deg3();

	test_degree_elevation_circle_deg2();
	test_degree_elevation_from_deg1();
	test_degree_elevation_from_deg2();

	test_is_closed();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////