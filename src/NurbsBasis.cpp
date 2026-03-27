#include "NurbsBasis.h"
#include "NurbsConstants.h"

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
            if (std::fabs(ndu[j][r]) > NurbsConstants::EpsilonNumerical)
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
            if (std::fabs(ndu[j][r]) > NurbsConstants::EpsilonNumerical)
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
            if (std::fabs(den) > NurbsConstants::EpsilonNumerical)
                temp = N[r] / den;
            N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}

} // namespace NurbsBasis