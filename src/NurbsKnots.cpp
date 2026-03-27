#include "NurbsKnots.h"

#include <algorithm>
#include <cmath>

namespace NurbsKnots {

void normalize_to_01(std::vector<double>& knots)
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

} // namespace NurbsKnots