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

// Options controlling intersection computation
struct IntersectionOptions
{
	double geomTol = 1e-6; // absolute geometric tolerance
	double paramTol = 1e-8; // parameter tolerance
	int seedSamplingRes = 16; // coarse grid resolution for seed detection
	double seedTolScale = 4.0; // seed threshold multiplier of geomTol
	int maxNewtonIterations = 30;
	double initialDamping = 1e-3; // LM initial damping
	bool verbose = false;
	double initialStep = 0.02; // predictor step in parameter-space
};

// Compute intersection between two NURBS surfaces
void compute_surface_intersection(const NurbsSurface& surface1, const NurbsSurface& surface2,
								NurbsIntersectionResult& result,
								const IntersectionOptions& options = IntersectionOptions());

#endif
