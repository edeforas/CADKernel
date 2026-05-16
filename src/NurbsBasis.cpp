#include "NurbsBasis.h"

#include <cmath>
#include <algorithm>

namespace NurbsBasis {

void basis_function_derivatives(
    int degree, 
    const std::vector<double>& knots, 
    int span, 
    double u, 
    std::vector<std::vector<double>>& ders)
{
    ders.assign(3, std::vector<double>(degree + 1, 0.));
    if (degree < 0)
        return;

    std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1, 0.));
    std::vector<double> left(degree + 1, 0.);
    std::vector<double> right(degree + 1, 0.);
    ndu[0][0] = 1.;

    for (int j = 1; j <= degree; ++j)
    {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        double saved = 0.;
        for (int r = 0; r < j; ++r)
        {
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = 0.;
            if (std::fabs(ndu[j][r]) > EpsilonNumerical)
                temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= degree; ++j)
        ders[0][j] = ndu[j][degree];

    std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1, 0.));
    for (int r = 0; r <= degree; ++r)
    {
        int s1 = 0;
        int s2 = 1;
        a[0][0] = 1.;

        for (int k = 1; k <= 2; ++k)
        {
            double d = 0.;
            int rk = r - k;
            int pk = degree - k;

            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? (k - 1) : (degree - r);
            for (int j = j1; j <= j2; ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }

    int r = degree;
    for (int k = 1; k <= 2; ++k)
    {
        for (int j = 0; j <= degree; ++j)
            ders[k][j] *= r;
        r *= (degree - k);
    }
}

void basis_function_derivatives(
    int degree, 
    const std::vector<double>& knots, 
    int span, 
    double u, 
    double ders[3][32])
{
    // Initialize to zero
    for (int k = 0; k < 3; ++k)
        for (int j = 0; j < 32; ++j)
            ders[k][j] = 0.;

    if (degree < 0)
        return;

    std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1, 0.));
    std::vector<double> left(degree + 1, 0.);
    std::vector<double> right(degree + 1, 0.);
    ndu[0][0] = 1.;

    for (int j = 1; j <= degree; ++j)
    {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        double saved = 0.;
        for (int r = 0; r < j; ++r)
        {
            ndu[j][r] = right[r + 1] + left[j - r];
            double temp = 0.;
            if (std::fabs(ndu[j][r]) > EpsilonNumerical)
                temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= degree; ++j)
        ders[0][j] = ndu[j][degree];

    std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1, 0.));
    for (int r = 0; r <= degree; ++r)
    {
        int s1 = 0;
        int s2 = 1;
        a[0][0] = 1.;

        for (int k = 1; k <= 2; ++k)
        {
            double d = 0.;
            int rk = r - k;
            int pk = degree - k;

            if (r >= k)
            {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            int j1 = (rk >= -1) ? 1 : -rk;
            int j2 = (r - 1 <= pk) ? (k - 1) : (degree - r);
            for (int j = j1; j <= j2; ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if (r <= pk)
            {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }

    int r = degree;
    for (int k = 1; k <= 2; ++k)
    {
        for (int j = 0; j <= degree; ++j)
            ders[k][j] *= r;
        r *= (degree - k);
    }
}

void basis_functions(
    int span, 
    double u, 
    int degree, 
    const std::vector<double>& knots, 
    std::vector<double>& N)
{
    N.assign(degree + 1, 0.);
    std::vector<double> left(degree + 1, 0.);
    std::vector<double> right(degree + 1, 0.);
    N[0] = 1.;

    for (int j = 1; j <= degree; ++j)
    {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        double saved = 0.;
        for (int r = 0; r < j; ++r)
        {
            const double den = right[r + 1] + left[j - r];
            double temp = 0.;
            if (std::fabs(den) > EpsilonNumerical)
                temp = N[r] / den;
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}



void normalize_to_01_knots(std::vector<double>& knots)
{
	if (knots.empty())
		return;

	//compute min and max
	double dMin = knots[0];
	double dMax = knots[0];
	for (size_t i = 1; i < knots.size(); ++i)
	{
		dMin = std::min(dMin, knots[i]);
		dMax = std::max(dMax, knots[i]);
	}

	if (dMax == dMin)
	{
		std::fill(knots.begin(), knots.end(), 0.0);
		return;
	}

	//apply scale so that 0<= knots <= 1
	for (auto& knot : knots)
		knot = (knot - dMin) / (dMax - dMin);
}

std::vector<double> build_uniform_knots(int degree, int nbCtrlPoints)
{
	std::vector<double> knots;
	if (nbCtrlPoints <= 0)
		return knots;

	if (degree < 0)
		degree = 0;
	if (degree >= nbCtrlPoints)
		degree = nbCtrlPoints - 1;

	knots.reserve(nbCtrlPoints + degree + 1);

	for (int i = 0; i <= degree; ++i)
		knots.push_back(0.);

	const int interior = nbCtrlPoints - degree - 1;
	for (int i = 1; i <= interior; ++i)
		knots.push_back((double)i);

	for (int i = 0; i <= degree; ++i)
		knots.push_back((double)(interior + 1));

	if (!knots.empty())
	{
		double dMin = knots.front();
		double dMax = knots.back();
		if (dMax > dMin)
			for (auto& k : knots)
				k = (k - dMin) / (dMax - dMin);
		else
			for (auto& k : knots)
				k = 0.;
	}

	return knots;
}

std::vector<double> build_clamped_uniform_knots(int degree, int nbPoints)
{
	std::vector<double> knots;
	if (nbPoints <= 0)
		return knots;

	if (degree < 0)
		degree = 0;

	if (degree >= nbPoints)
		degree = nbPoints - 1;

	knots.reserve(nbPoints + degree + 1);

	for (int i = 0; i <= degree; ++i)
		knots.push_back(0.);

	int interior = nbPoints - degree - 1;
	for (int i = 1; i <= interior; ++i)
		knots.push_back((double)i);

	for (int i = 0; i <= degree; ++i)
		knots.push_back((double)(interior + 1));

	return knots;
}

std::vector<double> build_open_uniform_knots(int degree, int nbPoints)
{
	std::vector<double> knots;
	if (nbPoints <= 0)
		return knots;

	if (degree < 1)
		degree = 1;
	if (degree >= nbPoints)
		degree = nbPoints - 1;

	knots.reserve(nbPoints + degree + 1);
	for (int i = 0; i <= degree; ++i)
		knots.push_back(0.);

	const int interior = nbPoints - degree - 1;
	for (int i = 1; i <= interior; ++i)
		knots.push_back((double)i / (double)(interior + 1));

	for (int i = 0; i <= degree; ++i)
		knots.push_back(1.);

	return knots;
}

std::vector<double> build_segmented_quadratic_knots(int nbSegments)
{
	std::vector<double> knots;
	if (nbSegments <= 0)
		return knots;

	knots.reserve(2 * nbSegments + 4);
	knots.push_back(0.);
	knots.push_back(0.);
	knots.push_back(0.);

	for (int i = 1; i < nbSegments; ++i)
	{
		const double t = (double)i / (double)nbSegments;
		knots.push_back(t);
		knots.push_back(t);
	}

	knots.push_back(1.);
	knots.push_back(1.);
	knots.push_back(1.);

	return knots;
}


} // namespace NurbsBasis