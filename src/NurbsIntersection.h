#ifndef NurbsIntersection_
#define NurbsIntersection_

#include "Geometry.h"

#include <vector>

struct NurbsIntersectionSample
{
	Point3 point;
	double uA;
	double vA;
	double uB;
	double vB;
};

struct NurbsIntersectionCurve
{
	std::vector<NurbsIntersectionSample> samples;
	bool closed;
};

struct NurbsIntersectionResult
{
	std::vector<NurbsIntersectionCurve> curves;
	bool hasPartialOverlap;
};

#endif
