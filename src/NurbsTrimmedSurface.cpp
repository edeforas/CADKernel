#include "NurbsTrimmedSurface.h"

#include <algorithm>
#include <cmath>

namespace
{
	double uv_dist_sq(const NurbsUvPoint& a, const NurbsUvPoint& b)
	{
		double du = a.u - b.u;
		double dv = a.v - b.v;
		return du * du + dv * dv;
	}

	double signed_area(const std::vector<NurbsUvPoint>& poly)
	{
		if (poly.size() < 3)
			return 0.;

		double area2 = 0.;
		for (int i = 0; i < (int)poly.size(); ++i)
		{
			const NurbsUvPoint& p0 = poly[i];
			const NurbsUvPoint& p1 = poly[(i + 1) % poly.size()];
			area2 += p0.u * p1.v - p1.u * p0.v;
		}
		return area2 * 0.5;
	}

	double orient2d(const NurbsUvPoint& a, const NurbsUvPoint& b, const NurbsUvPoint& c)
	{
		return (b.u - a.u) * (c.v - a.v) - (b.v - a.v) * (c.u - a.u);
	}

	bool segments_intersect_strict(const NurbsUvPoint& a, const NurbsUvPoint& b, const NurbsUvPoint& c, const NurbsUvPoint& d)
	{
		const double o1 = orient2d(a, b, c);
		const double o2 = orient2d(a, b, d);
		const double o3 = orient2d(c, d, a);
		const double o4 = orient2d(c, d, b);

		const double eps = 1.e-14;
		if (std::fabs(o1) <= eps || std::fabs(o2) <= eps || std::fabs(o3) <= eps || std::fabs(o4) <= eps)
			return false;

		return ((o1 > 0.) != (o2 > 0.)) && ((o3 > 0.) != (o4 > 0.));
	}

	bool has_self_intersection(const std::vector<NurbsUvPoint>& poly)
	{
		const int n = (int)poly.size();
		if (n < 4)
			return false;

		for (int i = 0; i < n; ++i)
		{
			const int i2 = (i + 1) % n;
			for (int j = i + 1; j < n; ++j)
			{
				const int j2 = (j + 1) % n;

				if (i == j || i2 == j || j2 == i)
					continue;

				if (segments_intersect_strict(poly[i], poly[i2], poly[j], poly[j2]))
					return true;
			}
		}

		return false;
	}

	void sort_by_angle_around_centroid(std::vector<NurbsUvPoint>& poly)
	{
		if (poly.size() < 3)
			return;

		double cu = 0.;
		double cv = 0.;
		for (const auto& p : poly)
		{
			cu += p.u;
			cv += p.v;
		}
		cu /= (double)poly.size();
		cv /= (double)poly.size();

		std::sort(poly.begin(), poly.end(), [cu, cv](const NurbsUvPoint& a, const NurbsUvPoint& b)
		{
			double aa = std::atan2(a.v - cv, a.u - cu);
			double ab = std::atan2(b.v - cv, b.u - cu);
			return aa < ab;
		});
	}

	void remove_collinear(std::vector<NurbsUvPoint>& poly)
	{
		if (poly.size() < 3)
			return;

		const double eps = 1.e-12;
		bool changed = true;
		while (changed && poly.size() >= 3)
		{
			changed = false;
			for (int i = 0; i < (int)poly.size(); ++i)
			{
				int i0 = (i - 1 + (int)poly.size()) % (int)poly.size();
				int i1 = i;
				int i2 = (i + 1) % (int)poly.size();

				if (uv_dist_sq(poly[i0], poly[i1]) <= eps || uv_dist_sq(poly[i1], poly[i2]) <= eps)
				{
					poly.erase(poly.begin() + i1);
					changed = true;
					break;
				}

				double cross = std::fabs(orient2d(poly[i0], poly[i1], poly[i2]));
				if (cross <= eps)
				{
					poly.erase(poly.begin() + i1);
					changed = true;
					break;
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////
NurbsTrimmedSurface::NurbsTrimmedSurface() :NurbsSurface()
{
}
/////////////////////////////////////////////////////////////////////////
NurbsTrimmedSurface::NurbsTrimmedSurface(const NurbsSurface& n) :NurbsSurface(n)
{
}
/////////////////////////////////////////////////////////////////////////
NurbsTrimmedSurface::NurbsTrimmedSurface(const NurbsTrimmedSurface& n) :
	NurbsSurface(n),
	_trimLoops(n._trimLoops)
{
}
/////////////////////////////////////////////////////////////////////////
NurbsTrimmedSurface::~NurbsTrimmedSurface()
{
}
/////////////////////////////////////////////////////////////////////////
NurbsTrimmedSurface& NurbsTrimmedSurface::operator=(const NurbsTrimmedSurface& other)
{
	if (this == &other)
		return *this;

	NurbsSurface::operator=(other);
	_trimLoops = other._trimLoops;
	return *this;
}
/////////////////////////////////////////////////////////////////////////
void NurbsTrimmedSurface::add_outer_loop(const std::vector<NurbsUvPoint>& loopPoints)
{
	add_loop(loopPoints, false);
}
/////////////////////////////////////////////////////////////////////////
void NurbsTrimmedSurface::add_inner_loop(const std::vector<NurbsUvPoint>& loopPoints)
{
	add_loop(loopPoints, true);
}
/////////////////////////////////////////////////////////////////////////
void NurbsTrimmedSurface::add_loop(const std::vector<NurbsUvPoint>& loopPoints, bool bHole)
{
	std::vector<NurbsUvPoint> fixed = loopPoints;
	if (!normalize_loop(fixed, bHole))
		return;

	NurbsTrimLoop loop;
	loop.points = fixed;
	loop.hole = bHole;
	_trimLoops.push_back(loop);
}
/////////////////////////////////////////////////////////////////////////
bool NurbsTrimmedSurface::normalize_loop(std::vector<NurbsUvPoint>& loopPoints, bool bHole)
{
	const double eps = 1.e-12;

	if (loopPoints.size() < 3)
		return false;

	for (auto& p : loopPoints)
	{
		p.u = std::clamp(p.u, 0., 1.);
		p.v = std::clamp(p.v, 0., 1.);
	}

	std::vector<NurbsUvPoint> compact;
	compact.reserve(loopPoints.size());
	for (int i = 0; i < (int)loopPoints.size(); ++i)
	{
		if (!compact.empty() && uv_dist_sq(compact.back(), loopPoints[i]) <= eps)
			continue;
		compact.push_back(loopPoints[i]);
	}

	if (compact.size() >= 2 && uv_dist_sq(compact.front(), compact.back()) <= eps)
		compact.pop_back();

	if (compact.size() < 3)
		return false;

	remove_collinear(compact);
	if (compact.size() < 3)
		return false;

	if (has_self_intersection(compact))
	{
		sort_by_angle_around_centroid(compact);
		remove_collinear(compact);
		if (compact.size() < 3 || has_self_intersection(compact))
			return false;
	}

	double area = signed_area(compact);
	if (std::fabs(area) <= 1.e-14)
		return false;

	if (!bHole && area < 0.)
		std::reverse(compact.begin(), compact.end());
	if (bHole && area > 0.)
		std::reverse(compact.begin(), compact.end());

	loopPoints.swap(compact);
	return true;
}
/////////////////////////////////////////////////////////////////////////
const std::vector<NurbsTrimLoop>& NurbsTrimmedSurface::trim_loops() const
{
	return _trimLoops;
}
/////////////////////////////////////////////////////////////////////////
bool NurbsTrimmedSurface::is_trimmed() const
{
	return !_trimLoops.empty();
}
/////////////////////////////////////////////////////////////////////////