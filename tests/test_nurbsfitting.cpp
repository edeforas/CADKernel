#include "NurbsFitting.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "StepFile.h"
#include "NurbsUtil.h"
#include "OBJFile.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

bool file_exists(const string& filename)
{
	ifstream in(filename.c_str());
	return in.is_open();
}

void test_nurbs_fitting_surface_grid_least_squares()
{
	cout << endl << "test_nurbs_fitting_surface_grid_least_squares" << endl;

	const int nUSamples = 17;
	const int nVSamples = 13;
	vector<Point3> samples(nUSamples * nVSamples);

	for (int v = 0; v < nVSamples; ++v)
	{
		const double tv = (double)v / (double)(nVSamples - 1);
		const double y = 6. * (tv - 0.5);
		for (int u = 0; u < nUSamples; ++u)
		{
			const double tu = (double)u / (double)(nUSamples - 1);
			const double x = 6. * (tu - 0.5);
			const double z = 0.08 * x * x + 0.03 * y * y + 0.25 * std::sin(0.7 * x) * std::cos(0.6 * y);
			samples[v * nUSamples + u] = Point3(x, y, z);
		}
	}

	NurbsSurface fit;
	const bool ok = NurbsFitting::fit_surface_least_squares(samples, nUSamples, nVSamples, 3, 3, 7, 6, fit);
	test_bool(ok, "surface fitting should succeed");
	test_bool(fit.nb_points_u() == 7, "control points U should match request");
	test_bool(fit.nb_points_v() == 6, "control points V should match request");

	double mse = 0.;
	for (int v = 0; v < nVSamples; ++v)
	{
		const double tv = (double)v / (double)(nVSamples - 1);
		for (int u = 0; u < nUSamples; ++u)
		{
			const double tu = (double)u / (double)(nUSamples - 1);
			Point3 p;
			fit.evaluate(tu, tv, p);
			const Point3 d = p - samples[v * nUSamples + u];
			mse += d.norm_square();
		}
	}
	mse /= (double)(nUSamples * nVSamples);
	const double rmse = std::sqrt(mse);
	test_bool(rmse < 0.5, "fitted surface RMSE should be bounded");

	StepWriter sw;
	sw.open("test_nurbs_fitting_surface_grid_least_squares.step");
	test_bool(sw.is_open(), "step writer should open fitted surface output");
	sw.write(fit);
	sw.close();

	Mesh m;
	NurbsUtil::to_mesh(fit, m, 10);

	OBJWriter ow;
	ow.open("test_nurbs_fitting_surface_grid_least_squares.obj");
	ow.write(m);


	test_bool(file_exists("test_nurbs_fitting_surface_grid_least_squares.step"), "fitted surface step file should exist");
}

void test_nurbs_fitting_curve_least_squares()
{
	cout << endl << "test_nurbs_fitting_curve_least_squares" << endl;

	const int nSamples = 61;
	vector<Point3> samples(nSamples);
	for (int i = 0; i < nSamples; ++i)
	{
		const double t = (double)i / (double)(nSamples - 1);
		const double x = 10. * (t - 0.5);
		const double y = 0.5 * std::sin(1.2 * x) + 0.05 * x * x;
		samples[i] = Point3(x, y, 0.);
	}

	NurbsCurve fit;
	const bool ok = NurbsFitting::fit_curve_least_squares(samples, 3, 12, fit);
	test_bool(ok, "curve fitting should succeed");
	test_bool(fit.nb_points() == 12, "curve control points should match request");

	double mse = 0.;
	for (int i = 0; i < nSamples; ++i)
	{
		const double t = (double)i / (double)(nSamples - 1);
		Point3 p;
		fit.evaluate(t, p);
		const Point3 d = p - samples[i];
		mse += d.norm_square();
	}
	mse /= (double)nSamples;
	const double rmse = std::sqrt(mse);
	test_bool(rmse < 0.25, "fitted curve RMSE should be bounded");

	Point3 p0, p1;
	fit.evaluate(0., p0);
	fit.evaluate(1., p1);
	test_bool((p0 - samples.front()).norm() < 0.5, "curve start should stay near first sample");
	test_bool((p1 - samples.back()).norm() < 0.5, "curve end should stay near last sample");
}

void test_nurbs_fitting_curve_line_accuracy()
{
	cout << endl << "test_nurbs_fitting_curve_line_accuracy" << endl;

	const int nSamples = 31;
	vector<Point3> samples(nSamples);
	for (int i = 0; i < nSamples; ++i)
	{
		const double t = (double)i / (double)(nSamples - 1);
		samples[i] = Point3(4. * t - 2., -1. + 3. * t, 2. * t);
	}

	NurbsCurve fit;
	const bool ok = NurbsFitting::fit_curve_least_squares(samples, 1, 2, fit);
	test_bool(ok, "line fitting should succeed");
	test_bool(fit.degree() == 1, "line fitting degree should be 1");
	test_bool(fit.nb_points() == 2, "line fitting should have 2 control points");

	double mse = 0.;
	for (int i = 0; i < nSamples; ++i)
	{
		const double t = (double)i / (double)(nSamples - 1);
		Point3 p;
		fit.evaluate(t, p);
		mse += (p - samples[i]).norm_square();
	}
	mse /= (double)nSamples;
	test_bool(std::sqrt(mse) < 1.e-6, "line fitting should be nearly exact");
}

void test_nurbs_fitting_surface_plane_accuracy()
{
	cout << endl << "test_nurbs_fitting_surface_plane_accuracy" << endl;

	const int nUSamples = 11;
	const int nVSamples = 9;
	vector<Point3> samples(nUSamples * nVSamples);

	for (int v = 0; v < nVSamples; ++v)
	{
		const double tv = (double)v / (double)(nVSamples - 1);
		for (int u = 0; u < nUSamples; ++u)
		{
			const double tu = (double)u / (double)(nUSamples - 1);
			const double x = 2. * tu - 1.;
			const double y = 4. * tv - 2.;
			const double z = 0.5 * x - 0.25 * y + 1.;
			samples[v * nUSamples + u] = Point3(x, y, z);
		}
	}

	NurbsSurface fit;
	const bool ok = NurbsFitting::fit_surface_least_squares(samples, nUSamples, nVSamples, 1, 1, 2, 2, fit);
	test_bool(ok, "plane fitting should succeed");
	test_bool(fit.degree_u() == 1 && fit.degree_v() == 1, "plane fitting degrees should be linear");

	double mse = 0.;
	for (int v = 0; v < nVSamples; ++v)
		for (int u = 0; u < nUSamples; ++u)
		{
			Point3 p;
			fit.evaluate((double)u / (double)(nUSamples - 1), (double)v / (double)(nVSamples - 1), p);
			mse += (p - samples[v * nUSamples + u]).norm_square();
		}
	mse /= (double)(nUSamples * nVSamples);
	test_bool(std::sqrt(mse) < 1.e-6, "plane fitting should be nearly exact");
}

void test_nurbs_fitting_invalid_inputs()
{
	cout << endl << "test_nurbs_fitting_invalid_inputs" << endl;

	NurbsCurve c;
	NurbsSurface s;

	test_bool(!NurbsFitting::fit_curve_least_squares({}, 3, 4, c), "curve fit should fail on empty samples");
	test_bool(!NurbsFitting::fit_curve_least_squares({ Point3(0.,0.,0.) }, 3, 2, c), "curve fit should fail on too few samples");
	test_bool(!NurbsFitting::fit_curve_least_squares({ Point3(0.,0.,0.), Point3(1.,0.,0.) }, 3, 3, c), "curve fit should fail when ctrl > samples");

	std::vector<Point3> grid(12, Point3(0., 0., 0.));
	test_bool(!NurbsFitting::fit_surface_least_squares(grid, 3, 4, 3, 3, 4, 4, s), "surface fit should fail when ctrl > samples");
	test_bool(!NurbsFitting::fit_surface_least_squares(grid, 3, 3, 2, 2, 2, 2, s), "surface fit should fail on mismatched sample size");
}

int main()
{
	test_nurbs_fitting_curve_least_squares();
	test_nurbs_fitting_curve_line_accuracy();
	test_nurbs_fitting_surface_grid_least_squares();
	test_nurbs_fitting_surface_plane_accuracy();
	test_nurbs_fitting_invalid_inputs();

	cout << "Test Finished.";
	return 0;
}
