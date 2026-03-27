#ifndef NurbsBoolean_
#define NurbsBoolean_

class NurbsSolid;
class NurbsTrimmedSurface;
struct NurbsIntersectionResult;

#include <vector>

class NurbsBoolean
{
public:
	NurbsBoolean();

	bool compute_union(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics = 0) const;
	bool compute_intersection(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics = 0) const;
	bool compute_difference(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics = 0) const;

	bool boolean_union(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const;
	bool boolean_intersection(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const;
	bool boolean_difference_bbox(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const;

};

#endif
