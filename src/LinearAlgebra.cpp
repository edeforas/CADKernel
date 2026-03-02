
#include "LinearAlgebra.h"
#include <cmath>

bool solve_linear_system(std::vector<double> M, const std::vector<double>& b, std::vector<double>& x, int n)
{
	x = b;
	if (n <= 0 || (int)M.size() != n * n || (int)b.size() != n)
		return false;

	for (int k = 0; k < n; ++k)
	{
		int iPivot = k;
		double dMax = std::fabs(M[k * n + k]);
		for (int i = k + 1; i < n; ++i)
		{
			const double v = std::fabs(M[i * n + k]);
			if (v > dMax)
			{
				dMax = v;
				iPivot = i;
			}
		}

		if (dMax < 1.e-20)
			return false;

		if (iPivot != k)
		{
			for (int j = k; j < n; ++j)
				std::swap(M[k * n + j], M[iPivot * n + j]);
			std::swap(x[k], x[iPivot]);
		}

		const double diag = M[k * n + k];
		for (int i = k + 1; i < n; ++i)
		{
			const double f = M[i * n + k] / diag;
			M[i * n + k] = 0.;
			for (int j = k + 1; j < n; ++j)
				M[i * n + j] -= f * M[k * n + j];
			x[i] -= f * x[k];
		}
	}

	for (int i = n - 1; i >= 0; --i)
	{
		double sum = x[i];
		for (int j = i + 1; j < n; ++j)
			sum -= M[i * n + j] * x[j];
		const double diag = M[i * n + i];
		if (std::fabs(diag) < 1.e-20)
			return false;
		x[i] = sum / diag;
	}

	return true;
}