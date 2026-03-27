#include "NurbsFitting.h"

#include "Geometry.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "LinearAlgebra.h"
#include "NurbsBasis.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace
{

int find_span(int nCtrl, int degree, double u, const std::vector<double>& knots)
{
	if (nCtrl <= 0 || knots.empty())
		return 0;

	if (u <= knots[degree])
		return degree;
	if (u >= knots[nCtrl])
		return nCtrl - 1;

	int low = degree;
	int high = nCtrl;
	int mid = (low + high) / 2;
	while (u < knots[mid] || u >= knots[mid + 1])
	{
		if (u < knots[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}

	return mid;
}

std::vector<double> build_basis_matrix(int nSamples, int nCtrl, int degree, const std::vector<double>& knots)
{
	std::vector<double> A(std::max(0, nSamples * nCtrl), 0.);
	if (nSamples <= 0 || nCtrl <= 0)
		return A;

	std::vector<double> N;
	for (int s = 0; s < nSamples; ++s)
	{
		double u = 0.;
		if (nSamples > 1)
			u = (double)s / (double)(nSamples - 1);

		const int span = find_span(nCtrl, degree, u, knots);
		NurbsBasis::basis_functions(span, u, degree, knots, N);
		const int first = span - degree;
		for (int j = 0; j <= degree; ++j)
		{
			const int iCtrl = first + j;
			if (iCtrl >= 0 && iCtrl < nCtrl)
				A[s * nCtrl + iCtrl] = N[j];
		}
	}

	return A;
}

std::vector<double> compute_AtA(const std::vector<double>& A, int m, int n)
{
	std::vector<double> AtA(n * n, 0.);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			double sum = 0.;
			for (int r = 0; r < m; ++r)
				sum += A[r * n + i] * A[r * n + j];
			AtA[i * n + j] = sum;
		}
	return AtA;
}

}

bool NurbsFitting::fit_curve_least_squares(
	const std::vector<Point3>& samples,
	int degree,
	int iNbCtrl,
	NurbsCurve& outCurve,
	double dRegularization)
{
	outCurve.clear();

	const int iNbSamples = (int)samples.size();
	if (iNbSamples < 2)
		return false;
	if (iNbCtrl < 2 || iNbCtrl > iNbSamples)
		return false;
	if (!std::isfinite(dRegularization) || dRegularization < 0.)
		return false;

	degree = std::max(1, std::min(degree, iNbCtrl - 1));

	const std::vector<double> knots = NurbsCurveUtil::build_open_uniform_knots(degree, iNbCtrl);
	if (knots.empty())
		return false;

	const std::vector<double> A = build_basis_matrix(iNbSamples, iNbCtrl, degree, knots);
	if (A.empty())
		return false;

	std::vector<double> M = compute_AtA(A, iNbSamples, iNbCtrl);
	for (int i = 0; i < iNbCtrl; ++i)
		M[i * iNbCtrl + i] += dRegularization;

	std::vector<double> bx(iNbCtrl, 0.);
	std::vector<double> by(iNbCtrl, 0.);
	std::vector<double> bz(iNbCtrl, 0.);

	for (int i = 0; i < iNbCtrl; ++i)
	{
		double sx = 0.;
		double sy = 0.;
		double sz = 0.;
		for (int s = 0; s < iNbSamples; ++s)
		{
			const double a = A[s * iNbCtrl + i];
			sx += a * samples[s].x();
			sy += a * samples[s].y();
			sz += a * samples[s].z();
		}
		bx[i] = sx;
		by[i] = sy;
		bz[i] = sz;
	}

	std::vector<double> sol;
	if (!solve_linear_system(M, bx, sol, iNbCtrl))
		return false;
	std::vector<Point3> ctrl(iNbCtrl, Point3(0., 0., 0.));
	for (int i = 0; i < iNbCtrl; ++i)
		ctrl[i].x() = sol[i];

	if (!solve_linear_system(M, by, sol, iNbCtrl))
		return false;
	for (int i = 0; i < iNbCtrl; ++i)
		ctrl[i].y() = sol[i];

	if (!solve_linear_system(M, bz, sol, iNbCtrl))
		return false;
	for (int i = 0; i < iNbCtrl; ++i)
		ctrl[i].z() = sol[i];

	outCurve.set_degree(degree);
	outCurve.set_points(ctrl);
	outCurve.set_knots(knots);
	outCurve.set_equals_weights();
	return true;
}

bool NurbsFitting::fit_surface_least_squares(
	const std::vector<Point3>& samples,
	int iNbSamplesU,
	int iNbSamplesV,
	int degreeU,
	int degreeV,
	int iNbCtrlU,
	int iNbCtrlV,
	NurbsSurface& outSurface,
	double dRegularization)
{
	outSurface.clear();

	if (iNbSamplesU < 2 || iNbSamplesV < 2)
		return false;
	if ((int)samples.size() != iNbSamplesU * iNbSamplesV)
		return false;
	if (iNbCtrlU < 2 || iNbCtrlV < 2)
		return false;
	if (iNbCtrlU > iNbSamplesU || iNbCtrlV > iNbSamplesV)
		return false;
	if (!std::isfinite(dRegularization) || dRegularization < 0.)
		return false;

	degreeU = std::max(1, std::min(degreeU, iNbCtrlU - 1));
	degreeV = std::max(1, std::min(degreeV, iNbCtrlV - 1));

	const std::vector<double> knotsU = NurbsCurveUtil::build_open_uniform_knots(degreeU, iNbCtrlU);
	const std::vector<double> knotsV = NurbsCurveUtil::build_open_uniform_knots(degreeV, iNbCtrlV);
	if (knotsU.empty() || knotsV.empty())
		return false;

	const std::vector<double> Au = build_basis_matrix(iNbSamplesU, iNbCtrlU, degreeU, knotsU);
	const std::vector<double> Av = build_basis_matrix(iNbSamplesV, iNbCtrlV, degreeV, knotsV);
	if (Au.empty() || Av.empty())
		return false;

	std::vector<double> Mu = compute_AtA(Au, iNbSamplesU, iNbCtrlU);
	std::vector<double> Mv = compute_AtA(Av, iNbSamplesV, iNbCtrlV);

	for (int i = 0; i < iNbCtrlU; ++i)
		Mu[i * iNbCtrlU + i] += dRegularization;
	for (int i = 0; i < iNbCtrlV; ++i)
		Mv[i * iNbCtrlV + i] += dRegularization;

	std::vector<double> bx(iNbCtrlU * iNbCtrlV, 0.);
	std::vector<double> by(iNbCtrlU * iNbCtrlV, 0.);
	std::vector<double> bz(iNbCtrlU * iNbCtrlV, 0.);

	for (int i = 0; i < iNbCtrlU; ++i)
		for (int j = 0; j < iNbCtrlV; ++j)
		{
			double sx = 0.;
			double sy = 0.;
			double sz = 0.;
			for (int su = 0; su < iNbSamplesU; ++su)
				for (int sv = 0; sv < iNbSamplesV; ++sv)
				{
					const double a = Au[su * iNbCtrlU + i] * Av[sv * iNbCtrlV + j];
					const Point3& p = samples[sv * iNbSamplesU + su];
					sx += a * p.x();
					sy += a * p.y();
					sz += a * p.z();
				}
			bx[j * iNbCtrlU + i] = sx;
			by[j * iNbCtrlU + i] = sy;
			bz[j * iNbCtrlU + i] = sz;
		}

	std::vector<double> tmpx(iNbCtrlU * iNbCtrlV, 0.);
	std::vector<double> tmpy(iNbCtrlU * iNbCtrlV, 0.);
	std::vector<double> tmpz(iNbCtrlU * iNbCtrlV, 0.);
	std::vector<double> rhs;
	std::vector<double> sol;

	for (int j = 0; j < iNbCtrlV; ++j)
	{
		rhs.assign(iNbCtrlU, 0.);
		for (int i = 0; i < iNbCtrlU; ++i) rhs[i] = bx[j * iNbCtrlU + i];
		if (!solve_linear_system(Mu, rhs, sol, iNbCtrlU)) return false;
		for (int i = 0; i < iNbCtrlU; ++i) tmpx[j * iNbCtrlU + i] = sol[i];

		for (int i = 0; i < iNbCtrlU; ++i) rhs[i] = by[j * iNbCtrlU + i];
		if (!solve_linear_system(Mu, rhs, sol, iNbCtrlU)) return false;
		for (int i = 0; i < iNbCtrlU; ++i) tmpy[j * iNbCtrlU + i] = sol[i];

		for (int i = 0; i < iNbCtrlU; ++i) rhs[i] = bz[j * iNbCtrlU + i];
		if (!solve_linear_system(Mu, rhs, sol, iNbCtrlU)) return false;
		for (int i = 0; i < iNbCtrlU; ++i) tmpz[j * iNbCtrlU + i] = sol[i];
	}

	std::vector<Point3> ctrl(iNbCtrlU * iNbCtrlV, Point3(0., 0., 0.));

	for (int i = 0; i < iNbCtrlU; ++i)
	{
		rhs.assign(iNbCtrlV, 0.);
		for (int j = 0; j < iNbCtrlV; ++j) rhs[j] = tmpx[j * iNbCtrlU + i];
		if (!solve_linear_system(Mv, rhs, sol, iNbCtrlV)) return false;
		for (int j = 0; j < iNbCtrlV; ++j) ctrl[i * iNbCtrlV + j].x() = sol[j];

		for (int j = 0; j < iNbCtrlV; ++j) rhs[j] = tmpy[j * iNbCtrlU + i];
		if (!solve_linear_system(Mv, rhs, sol, iNbCtrlV)) return false;
		for (int j = 0; j < iNbCtrlV; ++j) ctrl[i * iNbCtrlV + j].y() = sol[j];

		for (int j = 0; j < iNbCtrlV; ++j) rhs[j] = tmpz[j * iNbCtrlU + i];
		if (!solve_linear_system(Mv, rhs, sol, iNbCtrlV)) return false;
		for (int j = 0; j < iNbCtrlV; ++j) ctrl[i * iNbCtrlV + j].z() = sol[j];
	}

	std::vector<double> weights(iNbCtrlU * iNbCtrlV, 1.);
	outSurface.set_degree(degreeU, degreeV);
	outSurface.set_points(ctrl, iNbCtrlU, iNbCtrlV);
	outSurface.set_knots_u(knotsU);
	outSurface.set_knots_v(knotsV);
	outSurface.set_weights(weights);
	outSurface.set_closed_u(false);
	outSurface.set_closed_v(false);
	return true;
}
