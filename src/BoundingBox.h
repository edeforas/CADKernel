#ifndef BoundingBox_
#define BoundingBox_
#include "Geometry.h"

class BoundingBox
{
public:
	double xMin, yMin, zMin;
	double xMax, yMax, zMax;
	bool initialized;

	BoundingBox();

	void add(const Point3& p);

	bool intersects(const BoundingBox& b) const;
	bool contains(const BoundingBox& b) const;

	Point3 center() const;
	double diagonal_length() const;
};

enum class BBoxRelation
{
	EmptyA,
	EmptyB,
	Disjoint,
	AContainsB,
	BContainsA,
	PartialOverlap
};

double bbox_relation_tolerance(const BoundingBox& a, const BoundingBox& b);

bool bbox_intersects_tol(const BoundingBox& a, const BoundingBox& b, double tol);

bool bbox_contains_tol(const BoundingBox& outer, const BoundingBox& inner, double tol);

BBoxRelation classify_bbox_relation(const BoundingBox& a, const BoundingBox& b);

#endif