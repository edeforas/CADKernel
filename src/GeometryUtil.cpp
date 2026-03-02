#include "GeometryUtil.h"

#include <cmath>
#include <random>

namespace GeometryUtil
{
	void compute_random_axis(Point3& p, int iSeed)
	{
		std::default_random_engine generator;
		generator.seed(iSeed);
		std::normal_distribution<double> d; //use normal distribution to have points equally distributed on sphere  (instead of uniform)
		p.x() = d(generator);
		p.y() = d(generator);
		p.z() = d(generator);

		p.normalize();
	}
}