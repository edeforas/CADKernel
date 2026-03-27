#ifndef NurbsIntersection_
#define NurbsIntersection_

#include "Geometry.h"
#include "NurbsSurface.h"

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

// Compute intersection between two NURBS surfaces
void compute_surface_intersection(const NurbsSurface& surface1, const NurbsSurface& surface2,
                                NurbsIntersectionResult& result);

#endif
