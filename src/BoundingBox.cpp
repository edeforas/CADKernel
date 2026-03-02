#include "BoundingBox.h"

#include <algorithm>
#include <cmath>

#include "Geometry.h"


	BoundingBox::BoundingBox()
     : xMin(0.), yMin(0.), zMin(0.), xMax(0.), yMax(0.), zMax(0.), initialized(false)
     {}
	

	void BoundingBox::add(const Point3& p)
	{
		if (!initialized)
		{
			xMin = xMax = p.x();
			yMin = yMax = p.y();
			zMin = zMax = p.z();
			initialized = true;
			return;
		}

		xMin = std::min(xMin, p.x());
		yMin = std::min(yMin, p.y());
		zMin = std::min(zMin, p.z());
		xMax = std::max(xMax, p.x());
		yMax = std::max(yMax, p.y());
		zMax = std::max(zMax, p.z());
	}

	bool BoundingBox::intersects(const BoundingBox& b) const
	{
		if (!initialized || !b.initialized)
			return false;

		if (xMax < b.xMin || b.xMax < xMin)
			return false;
		if (yMax < b.yMin || b.yMax < yMin)
			return false;
		if (zMax < b.zMin || b.zMax < zMin)
			return false;

		return true;
	}

	bool BoundingBox::contains(const BoundingBox& b) const
	{
		if (!initialized || !b.initialized)
			return false;

		return (xMin <= b.xMin) && (yMin <= b.yMin) && (zMin <= b.zMin)
			&& (xMax >= b.xMax) && (yMax >= b.yMax) && (zMax >= b.zMax);
	}

	Point3 BoundingBox::center() const
	{
		return Point3((xMin + xMax) / 2., (yMin + yMax) / 2., (zMin + zMax) / 2.);
	}

	double BoundingBox::diagonal_length() const
	{
		if (!initialized)
			return 0.;

		double dx = xMax - xMin;
		double dy = yMax - yMin;
		double dz = zMax - zMin;
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

double bbox_relation_tolerance(const BoundingBox& a, const BoundingBox& b)
{
	const double dScale = std::max(a.diagonal_length(), b.diagonal_length());
	if (dScale <= 0.)
		return 1.e-12;
	return dScale * 1.e-9;
}

bool bbox_intersects_tol(const BoundingBox& a, const BoundingBox& b, double tol)
{
	if (!a.initialized || !b.initialized)
		return false;

	if (a.xMax < b.xMin - tol || b.xMax < a.xMin - tol)
		return false;
	if (a.yMax < b.yMin - tol || b.yMax < a.yMin - tol)
		return false;
	if (a.zMax < b.zMin - tol || b.zMax < a.zMin - tol)
		return false;

	return true;
}

bool bbox_contains_tol(const BoundingBox& outer, const BoundingBox& inner, double tol)
{
	if (!outer.initialized || !inner.initialized)
		return false;

	return (outer.xMin <= inner.xMin + tol) && (outer.yMin <= inner.yMin + tol) && (outer.zMin <= inner.zMin + tol)
		&& (outer.xMax >= inner.xMax - tol) && (outer.yMax >= inner.yMax - tol) && (outer.zMax >= inner.zMax - tol);
}

BBoxRelation classify_bbox_relation(const BoundingBox& a, const BoundingBox& b)
{
	if (!a.initialized)
		return BBoxRelation::EmptyA;
	if (!b.initialized)
		return BBoxRelation::EmptyB;

	const double tol = bbox_relation_tolerance(a, b);

	if (bbox_intersects_tol(a, b, tol))
		return BBoxRelation::PartialOverlap;

	if (bbox_contains_tol(a, b, tol))
		return BBoxRelation::AContainsB;

	if (bbox_contains_tol(b, a, tol))
		return BBoxRelation::BContainsA;

	return BBoxRelation::Disjoint;
}