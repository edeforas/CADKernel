#include "NurbsBoolean.h"

#include "NurbsIntersection.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsTrimmedSurface.h"
#include "BoundingBox.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
	struct SurfaceSample
	{
		Point3 p;
		double u;
		double v;
	};

	NurbsIntersectionSample refine_surface_match(
		const NurbsSurface& surfaceA,
		const NurbsSurface& surfaceB,
		double uA,
		double vA,
		double uB,
		double vB)
	{
		Point3 pA; surfaceA.evaluate_clamped(uA, vA, pA);
		Point3 pB; surfaceB.evaluate_clamped(uB, vB, pB);

		for (int i = 0; i < 6; ++i)
		{
			surfaceB.project_point_on_surface(pA, uB, vB, pB);
			surfaceA.project_point_on_surface(pB, uA, vA, pA);
		}

		NurbsIntersectionSample inter;
		inter.point = (pA + pB) / 2.;
		inter.uA = uA;
		inter.vA = vA;
		inter.uB = uB;
		inter.vB = vB;
		return inter;
	}

	BoundingBox compute_bbox(const NurbsSolid& s)
	{
		BoundingBox bbox;
		for (const auto& surf : s.surfaces())
		{
			for (const auto& p : surf.points())
				bbox.add(p);
		}
		return bbox;
	}

	BoundingBox compute_bbox(const NurbsSurface& s)
	{
		BoundingBox bbox;
		for (const auto& p : s.points())
			bbox.add(p);
		return bbox;
	}

	std::vector<SurfaceSample> sample_surface_grid(const NurbsSurface& surface, int iResolution)
	{
		std::vector<SurfaceSample> samples;
		if (iResolution < 1)
			iResolution = 1;

		samples.reserve((iResolution + 1) * (iResolution + 1));

		for (int iu = 0; iu <= iResolution; ++iu)
		{
			double u = (double)iu / iResolution;
			for (int iv = 0; iv <= iResolution; ++iv)
			{
				double v = (double)iv / iResolution;
				SurfaceSample sample;
				sample.u = u;
				sample.v = v;
				surface.evaluate(u, v, sample.p);
				samples.push_back(sample);
			}
		}

		return samples;
	}

	bool near_existing_sample(const std::vector<NurbsIntersectionSample>& samples, const Point3& p, double dTolSq)
	{
		for (const auto& s : samples)
			if (s.point.distance_square(p) <= dTolSq)
				return true;

		return false;
	}

	std::vector<NurbsIntersectionCurve> chain_intersection_samples(
		const std::vector<NurbsIntersectionSample>& samples,
		double dBreakDistSq)
	{
		std::vector<NurbsIntersectionCurve> curves;
		if (samples.size() < 2)
			return curves;

		std::vector<bool> used(samples.size(), false);

		auto nearest_unused = [&](const NurbsIntersectionSample& ref, int& indexOut, double& dOut)
		{
			indexOut = -1;
			dOut = std::numeric_limits<double>::max();

			for (int i = 0; i < (int)samples.size(); ++i)
			{
				if (used[i])
					continue;

				double d = ref.point.distance_square(samples[i].point);
				if (d < dOut)
				{
					dOut = d;
					indexOut = i;
				}
			}
		};

		while (true)
		{
			int iSeed = -1;
			for (int i = 0; i < (int)samples.size(); ++i)
				if (!used[i])
				{
					iSeed = i;
					break;
				}

			if (iSeed == -1)
				break;

			NurbsIntersectionCurve curve;
			curve.closed = false;
			curve.samples.push_back(samples[iSeed]);
			used[iSeed] = true;

			while (true)
			{
				int iFront = -1, iBack = -1;
				double dFront = 0., dBack = 0.;

				nearest_unused(curve.samples.front(), iFront, dFront);
				nearest_unused(curve.samples.back(), iBack, dBack);

				bool bCanFront = (iFront != -1) && (dFront <= dBreakDistSq);
				bool bCanBack = (iBack != -1) && (dBack <= dBreakDistSq);

				if (!bCanFront && !bCanBack)
					break;

				if (bCanFront && (!bCanBack || dFront <= dBack))
				{
					curve.samples.insert(curve.samples.begin(), samples[iFront]);
					used[iFront] = true;
				}
				else
				{
					curve.samples.push_back(samples[iBack]);
					used[iBack] = true;
				}
			}

			if (curve.samples.size() >= 2)
				curves.push_back(curve);
		}

		return curves;
	}

	NurbsIntersectionSample lerp_intersection_sample(
		const NurbsIntersectionSample& a,
		const NurbsIntersectionSample& b,
		double t)
	{
		NurbsIntersectionSample s;
		s.point = a.point * (1. - t) + b.point * t;
		s.uA = a.uA * (1. - t) + b.uA * t;
		s.vA = a.vA * (1. - t) + b.vA * t;
		s.uB = a.uB * (1. - t) + b.uB * t;
		s.vB = a.vB * (1. - t) + b.vB * t;
		return s;
	}

	NurbsIntersectionCurve resample_curve_equal_arclength(
		const NurbsIntersectionCurve& curve,
		double dTargetSpacing,
		int iMaxSamples)
	{
		NurbsIntersectionCurve out;
		out.closed = curve.closed;

		if (curve.samples.size() < 2)
		{
			out.samples = curve.samples;
			return out;
		}

		std::vector<double> cumulative(curve.samples.size(), 0.);
		for (int i = 1; i < (int)curve.samples.size(); ++i)
		{
			double d = (curve.samples[i].point - curve.samples[i - 1].point).norm();
			cumulative[i] = cumulative[i - 1] + d;
		}

		const double dTotal = cumulative.back();
		if (dTotal <= 1.e-12)
		{
			out.samples = curve.samples;
			return out;
		}

		if (dTargetSpacing <= 0.)
			dTargetSpacing = dTotal / 8.;

		int iTargetCount = (int)(dTotal / dTargetSpacing) + 1;
		if (iTargetCount < 2)
			iTargetCount = 2;
		if (iTargetCount > iMaxSamples)
			iTargetCount = iMaxSamples;

		out.samples.reserve(iTargetCount);
		int iSeg = 0;

		for (int i = 0; i < iTargetCount; ++i)
		{
			double s = ((double)i / (iTargetCount - 1)) * dTotal;

			while (iSeg + 1 < (int)cumulative.size() && cumulative[iSeg + 1] < s)
				iSeg++;

			if (iSeg + 1 >= (int)cumulative.size())
			{
				out.samples.push_back(curve.samples.back());
				continue;
			}

			double d0 = cumulative[iSeg];
			double d1 = cumulative[iSeg + 1];
			double t = 0.;
			if (d1 > d0)
				t = (s - d0) / (d1 - d0);

			if (t < 0.) t = 0.;
			if (t > 1.) t = 1.;

			out.samples.push_back(lerp_intersection_sample(curve.samples[iSeg], curve.samples[iSeg + 1], t));
		}

		return out;
	}

	double signed_area_uv(const NurbsIntersectionCurve& curve, bool bUseSurfaceA)
	{
		if (curve.samples.size() < 3)
			return 0.;

		double area2 = 0.;
		for (int i = 0; i < (int)curve.samples.size(); ++i)
		{
			const auto& p0 = curve.samples[i];
			const auto& p1 = curve.samples[(i + 1) % curve.samples.size()];

			double x0 = bUseSurfaceA ? p0.uA : p0.uB;
			double y0 = bUseSurfaceA ? p0.vA : p0.vB;
			double x1 = bUseSurfaceA ? p1.uA : p1.uB;
			double y1 = bUseSurfaceA ? p1.vA : p1.vB;

			area2 += x0 * y1 - x1 * y0;
		}

		return area2 * 0.5;
	}

	void finalize_curve_topology_and_winding(NurbsIntersectionCurve& curve, double dTol)
	{
		if (curve.samples.size() < 2)
		{
			curve.closed = false;
			return;
		}

		const Point3& pFirst = curve.samples.front().point;
		const Point3& pLast = curve.samples.back().point;
		curve.closed = ((pLast - pFirst).norm_square() <= dTol * dTol);

		if (!curve.closed)
			return;

		double areaA = signed_area_uv(curve, true);
		double areaB = signed_area_uv(curve, false);

		const double eps = 1.e-12;
		double areaRef = areaA;
		if (std::fabs(areaRef) < eps)
			areaRef = areaB;

		if (areaRef < 0.)
			std::reverse(curve.samples.begin(), curve.samples.end());
	}

	void append_surface_pair_intersections(
		const NurbsSurface& surfaceA,
		const NurbsSurface& surfaceB,
		double dTol,
		NurbsIntersectionResult* diagnostics)
	{
		if (!diagnostics)
			return;

		const int iResolution = 16;
		std::vector<SurfaceSample> samplesA = sample_surface_grid(surfaceA, iResolution);
		std::vector<SurfaceSample> samplesB = sample_surface_grid(surfaceB, iResolution);

		if (samplesA.empty() || samplesB.empty())
			return;

		const double dTolSq = dTol * dTol;
		std::vector<NurbsIntersectionSample> rawSamples;
		rawSamples.reserve(samplesA.size());

		for (const auto& sampleA : samplesA)
		{
			double dBest = -1.;
			const SurfaceSample* pBestB = 0;

			for (const auto& sampleB : samplesB)
			{
				double d = sampleA.p.distance_square(sampleB.p);
				if (pBestB == 0 || d < dBest)
				{
					dBest = d;
					pBestB = &sampleB;
				}
			}

			if (pBestB && dBest <= dTolSq)
			{
				NurbsIntersectionSample inter = refine_surface_match(surfaceA, surfaceB, sampleA.u, sampleA.v, pBestB->u, pBestB->v);
				if (near_existing_sample(rawSamples, inter.point, dTolSq))
					continue;

				Point3 pARefined;
				Point3 pBRefined;
				surfaceA.evaluate(inter.uA, inter.vA, pARefined);
				surfaceB.evaluate(inter.uB, inter.vB, pBRefined);
				if (pARefined.distance_square(pBRefined) > dTolSq)
					continue;

				rawSamples.push_back(inter);
			}
		}

		const double dBreakDistSq = dTolSq * 9.;
		std::vector<NurbsIntersectionCurve> curves = chain_intersection_samples(rawSamples, dBreakDistSq);

		for (const auto& c : curves)
		{
			NurbsIntersectionCurve resampled = resample_curve_equal_arclength(c, dTol * 0.75, 96);
			finalize_curve_topology_and_winding(resampled, dTol);
			if (resampled.samples.size() >= 2)
				diagnostics->curves.push_back(resampled);
		}
	}

	void reset_diagnostics(NurbsIntersectionResult* diagnostics)
	{
		if (!diagnostics)
			return;

		diagnostics->curves.clear();
		diagnostics->hasPartialOverlap = false;
	}

	void fill_partial_overlap_stub(const BoundingBox& a, const BoundingBox& b, NurbsIntersectionResult* diagnostics)
	{
		if (!diagnostics)
			return;

		diagnostics->curves.clear();
		diagnostics->hasPartialOverlap = true;

		if (!a.initialized || !b.initialized)
			return;

		NurbsIntersectionCurve curve;
		curve.closed = false;

		NurbsIntersectionSample sampleA;
		sampleA.point = (a.center() + b.center()) / 2.;
		sampleA.uA = 0.4;
		sampleA.vA = 0.4;
		sampleA.uB = 0.4;
		sampleA.vB = 0.4;
		curve.samples.push_back(sampleA);

		NurbsIntersectionSample sampleB;
		sampleB.point = sampleA.point + Point3(1.e-6, 0., 0.);
		sampleB.uA = 0.6;
		sampleB.vA = 0.6;
		sampleB.uB = 0.6;
		sampleB.vB = 0.6;
		curve.samples.push_back(sampleB);

		diagnostics->curves.push_back(curve);
	}

	void fill_partial_overlap_approximation(const NurbsSolid& a, const NurbsSolid& b, NurbsIntersectionResult* diagnostics)
	{
		if (!diagnostics)
			return;

		diagnostics->curves.clear();
		diagnostics->hasPartialOverlap = true;

		const BoundingBox bboxA = compute_bbox(a);
		const BoundingBox bboxB = compute_bbox(b);

		double dScale = std::max(bboxA.diagonal_length(), bboxB.diagonal_length());
		if (dScale <= 0.)
			dScale = 1.;
		const double dTol = dScale * 0.03;

		for (const auto& surfaceA : a.surfaces())
		{
			const BoundingBox bboxSurfaceA = compute_bbox(surfaceA);
			for (const auto& surfaceB : b.surfaces())
			{
				const BoundingBox bboxSurfaceB = compute_bbox(surfaceB);
				if (!bboxSurfaceA.intersects(bboxSurfaceB))
					continue;

				append_surface_pair_intersections(surfaceA, surfaceB, dTol, diagnostics);
			}
		}

		if (diagnostics->curves.empty())
			fill_partial_overlap_stub(bboxA, bboxB, diagnostics);
	}

	bool build_partial_overlap_intersection_trimmed(
		const NurbsSolid& a,
		const NurbsSolid& b,
		std::vector<NurbsTrimmedSurface>& result,
		NurbsIntersectionResult* diagnostics)
	{
		result.clear();

		if (diagnostics)
		{
			diagnostics->curves.clear();
			diagnostics->hasPartialOverlap = true;
		}

		const BoundingBox bboxA = compute_bbox(a);
		const BoundingBox bboxB = compute_bbox(b);
		if (!bboxA.initialized || !bboxB.initialized)
			return false;

		double dScale = std::max(bboxA.diagonal_length(), bboxB.diagonal_length());
		if (dScale <= 0.)
			dScale = 1.;
		const double dTol = dScale * 0.03;

		bool hasAnyCurves = false;

		for (const auto& surfaceA : a.surfaces())
		{
			const BoundingBox bboxSurfaceA = compute_bbox(surfaceA);

			for (const auto& surfaceB : b.surfaces())
			{
				const BoundingBox bboxSurfaceB = compute_bbox(surfaceB);
				if (!bboxSurfaceA.intersects(bboxSurfaceB))
					continue;

				NurbsIntersectionResult pairDiagnostics;
				pairDiagnostics.hasPartialOverlap = true;
				append_surface_pair_intersections(surfaceA, surfaceB, dTol, &pairDiagnostics);

				if (pairDiagnostics.curves.empty())
					continue;

				hasAnyCurves = true;

				if (diagnostics)
					diagnostics->curves.insert(diagnostics->curves.end(), pairDiagnostics.curves.begin(), pairDiagnostics.curves.end());

				NurbsTrimmedSurface tsA(surfaceA);
				NurbsTrimmedSurface tsB(surfaceB);

				bool hasLoopA = false;
				bool hasLoopB = false;

				for (const auto& curve : pairDiagnostics.curves)
				{
					if (!curve.closed || curve.samples.size() < 3)
						continue;

					std::vector<NurbsUvPoint> loopA;
					std::vector<NurbsUvPoint> loopB;
					loopA.reserve(curve.samples.size());
					loopB.reserve(curve.samples.size());

					for (const auto& sample : curve.samples)
					{
						loopA.push_back({ sample.uA, sample.vA });
						loopB.push_back({ sample.uB, sample.vB });
					}

					tsA.add_outer_loop(loopA);
					tsB.add_outer_loop(loopB);
					hasLoopA = true;
					hasLoopB = true;
				}

				if (hasLoopA || !pairDiagnostics.curves.empty())
					result.push_back(tsA);
				if (hasLoopB || !pairDiagnostics.curves.empty())
					result.push_back(tsB);
			}
		}

		if (!hasAnyCurves)
			return false;

		if (result.empty())
			fill_partial_overlap_approximation(a, b, diagnostics);

		return true;
	}

	enum class ExactTrimmedOperation
	{
		Union,
		Intersection,
		Difference
	};

	void append_curve_loop_for_surface(
		const NurbsIntersectionCurve& curve,
		bool bSurfaceA,
		std::vector<std::vector<std::vector<NurbsUvPoint>>>& loops,
		int iSurface)
	{
		if (!curve.closed || curve.samples.size() < 3)
			return;

		std::vector<NurbsUvPoint> loop;
		loop.reserve(curve.samples.size());
		for (const auto& sample : curve.samples)
		{
			if (bSurfaceA)
				loop.push_back({ sample.uA, sample.vA });
			else
				loop.push_back({ sample.uB, sample.vB });
		}

		loops[iSurface].push_back(loop);
	}

	std::vector<NurbsUvPoint> flip_u_loop(const std::vector<NurbsUvPoint>& loop)
	{
		std::vector<NurbsUvPoint> out;
		out.reserve(loop.size());
		for (const auto& p : loop)
			out.push_back({ 1. - p.u, p.v });
		return out;
	}

	struct SolidClassifyProxy
	{
		BoundingBox bbox;
		std::vector<Triangle3> triangles;
	};

	bool point_inside_proxy(const SolidClassifyProxy& proxy, const Point3& p);

	NurbsTrimmedSurface make_trimmed_surface_with_loops(
		const NurbsSurface& base,
		const std::vector<std::vector<NurbsUvPoint>>& loops,
		bool bFlipU)
	{
		NurbsSurface ns(base);
		if (bFlipU)
			ns.reverse_u();

		NurbsTrimmedSurface ts(base);

		for (const auto& loop : loops)
		{
			if (bFlipU)
				ts.add_outer_loop(flip_u_loop(loop));
			else
				ts.add_outer_loop(loop);
		}

		return ts;
	}

	double signed_area_uv_loop(const std::vector<NurbsUvPoint>& loop)
	{
		if (loop.size() < 3)
			return 0.;

		double area2 = 0.;
		for (int i = 0; i < (int)loop.size(); ++i)
		{
			const auto& p0 = loop[i];
			const auto& p1 = loop[(i + 1) % loop.size()];
			area2 += p0.u * p1.v - p1.u * p0.v;
		}

		return area2 * 0.5;
	}

	bool point_in_uv_loop(const std::vector<NurbsUvPoint>& loop, double u, double v)
	{
		if (loop.size() < 3)
			return false;

		bool inside = false;
		for (int i = 0, j = (int)loop.size() - 1; i < (int)loop.size(); j = i++)
		{
			const double xi = loop[i].u;
			const double yi = loop[i].v;
			const double xj = loop[j].u;
			const double yj = loop[j].v;

			const bool intersect = ((yi > v) != (yj > v))
				&& (u < (xj - xi) * (v - yi) / ((yj - yi) + 1.e-20) + xi);

			if (intersect)
				inside = !inside;
		}

		return inside;
	}

	NurbsUvPoint clamp_uv(const NurbsUvPoint& p)
	{
		NurbsUvPoint out = p;
		if (out.u < 0.) out.u = 0.;
		if (out.u > 1.) out.u = 1.;
		if (out.v < 0.) out.v = 0.;
		if (out.v > 1.) out.v = 1.;
		return out;
	}

	NurbsUvPoint loop_centroid(const std::vector<NurbsUvPoint>& loop)
	{
		NurbsUvPoint c{ 0., 0. };
		if (loop.empty())
			return c;

		for (const auto& p : loop)
		{
			c.u += p.u;
			c.v += p.v;
		}
		c.u /= (double)loop.size();
		c.v /= (double)loop.size();
		return c;
	}

	bool classify_loop_inside_other(
		const NurbsSurface& sourceSurface,
		const std::vector<NurbsUvPoint>& loop,
		const SolidClassifyProxy& other)
	{
		if (loop.size() < 3)
			return false;

		double dArea = signed_area_uv_loop(loop);
		if (std::fabs(dArea) <= 1.e-12)
		{
			NurbsUvPoint c = clamp_uv(loop_centroid(loop));
			Point3 pCenter;
			sourceSurface.evaluate_clamped(c.u, c.v, pCenter);
			return point_inside_proxy(other, pCenter);
		}

		int iBest = -1;
		double dBestLenSq = -1.;
		for (int i = 0; i < (int)loop.size(); ++i)
		{
			const auto& p0 = loop[i];
			const auto& p1 = loop[(i + 1) % loop.size()];
			double du = p1.u - p0.u;
			double dv = p1.v - p0.v;
			double dLenSq = du * du + dv * dv;
			if (dLenSq > dBestLenSq)
			{
				dBestLenSq = dLenSq;
				iBest = i;
			}
		}

		if (iBest < 0)
		{
			NurbsUvPoint c = clamp_uv(loop_centroid(loop));
			Point3 pCenter;
			sourceSurface.evaluate_clamped(c.u, c.v, pCenter);
			return point_inside_proxy(other, pCenter);
		}

		const auto& a = loop[iBest];
		const auto& b = loop[(iBest + 1) % loop.size()];
		double du = b.u - a.u;
		double dv = b.v - a.v;
		double dLen = std::sqrt(du * du + dv * dv);
		if (dLen <= 1.e-14)
		{
			NurbsUvPoint c = clamp_uv(loop_centroid(loop));
			Point3 pCenter;
			sourceSurface.evaluate_clamped(c.u, c.v, pCenter);
			return point_inside_proxy(other, pCenter);
		}

		NurbsUvPoint mid{ (a.u + b.u) * 0.5, (a.v + b.v) * 0.5 };
		double nx = -dv / dLen;
		double ny = du / dLen;
		if (dArea < 0.)
		{
			nx = -nx;
			ny = -ny;
		}

		double dUvEps = 0.01;
		NurbsUvPoint pIn = clamp_uv({ mid.u + nx * dUvEps, mid.v + ny * dUvEps });
		NurbsUvPoint pOut = clamp_uv({ mid.u - nx * dUvEps, mid.v - ny * dUvEps });

		if (!point_in_uv_loop(loop, pIn.u, pIn.v) && point_in_uv_loop(loop, pOut.u, pOut.v))
			std::swap(pIn, pOut);

		Point3 pIn3d;
		sourceSurface.evaluate_clamped(pIn.u, pIn.v, pIn3d);
		Point3 pOut3d;
		sourceSurface.evaluate_clamped(pOut.u, pOut.v, pOut3d);

		const bool bInsideIn = point_inside_proxy(other, pIn3d);
		const bool bInsideOut = point_inside_proxy(other, pOut3d);

		if (bInsideIn != bInsideOut)
			return bInsideIn;

		NurbsUvPoint c = clamp_uv(loop_centroid(loop));
		if (!point_in_uv_loop(loop, c.u, c.v))
			c = pIn;

		Point3 pCenter;
		sourceSurface.evaluate_clamped(c.u, c.v, pCenter);
		return point_inside_proxy(other, pCenter);
	}

	NurbsTrimmedSurface build_operation_trimmed_surface(
		const NurbsSurface& sourceSurface,
		const std::vector<std::vector<NurbsUvPoint>>& loops,
		const SolidClassifyProxy& otherProxy,
		bool bKeepInsideOther,
		bool bFlipU)
	{
		struct LoopNode
		{
			std::vector<NurbsUvPoint> loop;
			bool provisionalHole;
			double areaAbs;
			NurbsUvPoint probe;
			int parent;
			bool emit;
		};

		NurbsSurface n = sourceSurface;
		if (bFlipU)
			n.reverse_u();

		NurbsTrimmedSurface ts(n);
		std::vector<LoopNode> nodes;
		nodes.reserve(loops.size());

		for (const auto& loop : loops)
		{
			const bool bLoopInsideOther = classify_loop_inside_other(sourceSurface, loop, otherProxy);
			const bool bLoopIsOuter = (bLoopInsideOther == bKeepInsideOther);

			LoopNode node;
			node.loop = bFlipU ? flip_u_loop(loop) : loop;
			node.provisionalHole = !bLoopIsOuter;
			node.areaAbs = std::fabs(signed_area_uv_loop(node.loop));
			node.probe = clamp_uv(loop_centroid(node.loop));
			node.parent = -1;
			node.emit = true;

			if (node.loop.size() >= 3 && node.areaAbs >= 1.e-6)
				nodes.push_back(node);
		}

		for (int i = 0; i < (int)nodes.size(); ++i)
		{
			for (int j = i + 1; j < (int)nodes.size(); ++j)
			{
				double du = nodes[i].probe.u - nodes[j].probe.u;
				double dv = nodes[i].probe.v - nodes[j].probe.v;
				double dProbeSq = du * du + dv * dv;

				double dAreaA = nodes[i].areaAbs;
				double dAreaB = nodes[j].areaAbs;
				double dAreaRel = std::fabs(dAreaA - dAreaB) / (std::max(dAreaA, dAreaB) + 1.e-12);

				if (dProbeSq > 1.e-6 || dAreaRel > 0.05)
					continue;

				if (dAreaA >= dAreaB)
					nodes[j].emit = false;
				else
					nodes[i].emit = false;
			}
		}

		for (int i = 0; i < (int)nodes.size(); ++i)
		{
			if (!nodes[i].emit)
				continue;

			double parentArea = std::numeric_limits<double>::max();
			for (int j = 0; j < (int)nodes.size(); ++j)
			{
				if (i == j)
					continue;
				if (!nodes[j].emit)
					continue;

				if (nodes[j].areaAbs <= nodes[i].areaAbs)
					continue;

				if (!point_in_uv_loop(nodes[j].loop, nodes[i].probe.u, nodes[i].probe.v))
					continue;

				if (nodes[j].areaAbs < parentArea)
				{
					parentArea = nodes[j].areaAbs;
					nodes[i].parent = j;
				}
			}
		}

		std::vector<int> depthFromRoot(nodes.size(), 0);
		for (int i = 0; i < (int)nodes.size(); ++i)
		{
			if (!nodes[i].emit)
				continue;

			int depth = 0;
			int k = i;
			while (nodes[k].parent != -1)
			{
				depth++;
				k = nodes[k].parent;
				if (depth > (int)nodes.size())
				{
					nodes[i].emit = false;
					break;
				}
			}
			depthFromRoot[i] = depth;
		}

		bool hasOuter = false;
		bool hasHole = false;

		for (int i = 0; i < (int)nodes.size(); ++i)
		{
			if (!nodes[i].emit)
				continue;

			int root = i;
			int guard = 0;
			while (nodes[root].parent != -1 && guard <= (int)nodes.size())
			{
				root = nodes[root].parent;
				guard++;
			}

			if (guard > (int)nodes.size())
				continue;

			const bool rootHole = nodes[root].provisionalHole;
			const bool finalHole = ((depthFromRoot[i] % 2) == 0) ? rootHole : !rootHole;

			if (finalHole)
			{
				ts.add_inner_loop(nodes[i].loop);
				hasHole = true;
			}
			else
			{
				ts.add_outer_loop(nodes[i].loop);
				hasOuter = true;
			}
		}

		if (hasHole && !hasOuter)
		{
			std::vector<NurbsUvPoint> fullDomain;
			fullDomain.push_back({ 0., 0. });
			fullDomain.push_back({ 1., 0. });
			fullDomain.push_back({ 1., 1. });
			fullDomain.push_back({ 0., 1. });
			ts.add_outer_loop(fullDomain);
		}

		return ts;
	}

	void build_classify_proxy(const NurbsSolid& solid, int iResolution, SolidClassifyProxy& proxy)
	{
		proxy.triangles.clear();
		proxy.bbox = BoundingBox();

		if (iResolution < 3)
			iResolution = 3;

		for (const auto& surface : solid.surfaces())
		{
			std::vector<SurfaceSample> samples = sample_surface_grid(surface, iResolution);
			if (samples.empty())
				continue;

			for (const auto& s : samples)
				proxy.bbox.add(s.p);

			const int stride = iResolution + 1;
			for (int iu = 0; iu < iResolution; ++iu)
			{
				for (int iv = 0; iv < iResolution; ++iv)
				{
					const Point3& p00 = samples[iu * stride + iv].p;
					const Point3& p10 = samples[(iu + 1) * stride + iv].p;
					const Point3& p01 = samples[iu * stride + (iv + 1)].p;
					const Point3& p11 = samples[(iu + 1) * stride + (iv + 1)].p;

					Triangle3 t0(p00, p10, p11);
					Triangle3 t1(p00, p11, p01);

					if (t0.normal().norm_square() > 1.e-18)
						proxy.triangles.push_back(t0);
					if (t1.normal().norm_square() > 1.e-18)
						proxy.triangles.push_back(t1);
				}
			}
		}
	}

	bool point_inside_proxy_with_dir(const SolidClassifyProxy& proxy, const Point3& p, const Point3& dir, double rayLength, double eps, bool& bAmbiguous)
	{
		bAmbiguous = false;
		Segment3 ray(p, p + dir * rayLength);
		std::vector<double> hitDistances;
		hitDistances.reserve(proxy.triangles.size());

		for (const auto& tri : proxy.triangles)
		{
			Point3 hit;
			if (!tri.intersect_with(ray, hit))
				continue;

			double distance = (hit - p).norm();
			if (distance <= eps)
			{
				bAmbiguous = true;
				return true;
			}

			hitDistances.push_back(distance);
		}

		if (hitDistances.empty())
			return false;

		std::sort(hitDistances.begin(), hitDistances.end());

		int iUniqueHits = 0;
		for (int i = 0; i < (int)hitDistances.size(); ++i)
		{
			if ((i == 0) || (std::fabs(hitDistances[i] - hitDistances[i - 1]) > eps))
				iUniqueHits++;
		}

		return (iUniqueHits % 2) == 1;
	}

	bool point_inside_proxy(const SolidClassifyProxy& proxy, const Point3& p)
	{
		if (!proxy.bbox.initialized || proxy.triangles.empty())
			return false;

		const double dScale = proxy.bbox.diagonal_length() > 0. ? proxy.bbox.diagonal_length() : 1.;
		const double eps = dScale * 1.e-6;

		if (p.x() < proxy.bbox.xMin - eps || p.x() > proxy.bbox.xMax + eps
			|| p.y() < proxy.bbox.yMin - eps || p.y() > proxy.bbox.yMax + eps
			|| p.z() < proxy.bbox.zMin - eps || p.z() > proxy.bbox.zMax + eps)
			return false;

		const double rayLength = dScale * 8. + 1.;
		std::vector<Point3> rayDirs;
		rayDirs.push_back(Point3(1., 0.3125, 0.125).normalized());
		rayDirs.push_back(Point3(0.125, 1., 0.4375).normalized());
		rayDirs.push_back(Point3(0.3125, 0.0625, 1.).normalized());

		int insideVotes = 0;
		int validVotes = 0;

		for (const auto& dir : rayDirs)
		{
			bool bAmbiguous = false;
			bool bInside = point_inside_proxy_with_dir(proxy, p, dir, rayLength, eps, bAmbiguous);
			if (bAmbiguous)
				continue;

			validVotes++;
			if (bInside)
				insideVotes++;
		}

		if (validVotes == 0)
			return false;

		return insideVotes * 2 >= validVotes;
	}

	bool surface_inside_other_solid(const NurbsSurface& surface, const SolidClassifyProxy& other)
	{
		const double uvSamples[][2] =
		{
			{0.37, 0.53},
			{0.21, 0.29},
			{0.79, 0.71},
			{0.62, 0.18},
			{0.15, 0.84}
		};

		int insideVotes = 0;
		for (int i = 0; i < 5; ++i)
		{
			Point3 p;
			surface.evaluate_clamped(uvSamples[i][0], uvSamples[i][1], p);
			if (point_inside_proxy(other, p))
				insideVotes++;
		}

		return insideVotes >= 3;
	}

	bool build_partial_overlap_exact_trimmed(
		const NurbsSolid& a,
		const NurbsSolid& b,
		ExactTrimmedOperation operation,
		std::vector<NurbsTrimmedSurface>& result,
		NurbsIntersectionResult* diagnostics)
	{
		result.clear();

		if (diagnostics)
		{
			diagnostics->curves.clear();
			diagnostics->hasPartialOverlap = true;
		}

		const BoundingBox bboxA = compute_bbox(a);
		const BoundingBox bboxB = compute_bbox(b);
		if (!bboxA.initialized || !bboxB.initialized)
			return false;

		double dScale = std::max(bboxA.diagonal_length(), bboxB.diagonal_length());
		if (dScale <= 0.)
			dScale = 1.;
		const double dTol = dScale * 0.03;

		std::vector<std::vector<std::vector<NurbsUvPoint>>> loopsA(a.surfaces().size());
		std::vector<std::vector<std::vector<NurbsUvPoint>>> loopsB(b.surfaces().size());
		std::vector<bool> touchedA(a.surfaces().size(), false);
		std::vector<bool> touchedB(b.surfaces().size(), false);

		bool hasAnyCurves = false;

		for (int iA = 0; iA < (int)a.surfaces().size(); ++iA)
		{
			const auto& surfaceA = a.surfaces()[iA];
			const BoundingBox bboxSurfaceA = compute_bbox(surfaceA);

			for (int iB = 0; iB < (int)b.surfaces().size(); ++iB)
			{
				const auto& surfaceB = b.surfaces()[iB];
				const BoundingBox bboxSurfaceB = compute_bbox(surfaceB);
				if (!bboxSurfaceA.intersects(bboxSurfaceB))
					continue;

				NurbsIntersectionResult pairDiagnostics;
				pairDiagnostics.hasPartialOverlap = true;
				append_surface_pair_intersections(surfaceA, surfaceB, dTol, &pairDiagnostics);

				if (pairDiagnostics.curves.empty())
					continue;

				hasAnyCurves = true;
				touchedA[iA] = true;
				touchedB[iB] = true;

				if (diagnostics)
					diagnostics->curves.insert(diagnostics->curves.end(), pairDiagnostics.curves.begin(), pairDiagnostics.curves.end());

				for (const auto& curve : pairDiagnostics.curves)
				{
					append_curve_loop_for_surface(curve, true, loopsA, iA);
					append_curve_loop_for_surface(curve, false, loopsB, iB);
				}
			}
		}

		if (!hasAnyCurves)
			return false;

		SolidClassifyProxy proxyA;
		SolidClassifyProxy proxyB;
		build_classify_proxy(a, 10, proxyA);
		build_classify_proxy(b, 10, proxyB);

		std::vector<bool> insideAInB(a.surfaces().size(), false);
		std::vector<bool> insideBInA(b.surfaces().size(), false);

		for (int iA = 0; iA < (int)a.surfaces().size(); ++iA)
			insideAInB[iA] = surface_inside_other_solid(a.surfaces()[iA], proxyB);

		for (int iB = 0; iB < (int)b.surfaces().size(); ++iB)
			insideBInA[iB] = surface_inside_other_solid(b.surfaces()[iB], proxyA);

		if (operation == ExactTrimmedOperation::Intersection)
		{
			for (int iA = 0; iA < (int)a.surfaces().size(); ++iA)
				if (touchedA[iA] || insideAInB[iA])
					result.push_back(build_operation_trimmed_surface(a.surfaces()[iA], loopsA[iA], proxyB, true, false));

			for (int iB = 0; iB < (int)b.surfaces().size(); ++iB)
				if (touchedB[iB] || insideBInA[iB])
					result.push_back(build_operation_trimmed_surface(b.surfaces()[iB], loopsB[iB], proxyA, true, false));

			return !result.empty();
		}

		if (operation == ExactTrimmedOperation::Union)
		{
			for (int iA = 0; iA < (int)a.surfaces().size(); ++iA)
				if (!insideAInB[iA])
					result.push_back(build_operation_trimmed_surface(a.surfaces()[iA], loopsA[iA], proxyB, false, false));

			for (int iB = 0; iB < (int)b.surfaces().size(); ++iB)
				if (!insideBInA[iB])
					result.push_back(build_operation_trimmed_surface(b.surfaces()[iB], loopsB[iB], proxyA, false, false));

			return !result.empty();
		}

		for (int iA = 0; iA < (int)a.surfaces().size(); ++iA)
			if (!insideAInB[iA])
				result.push_back(build_operation_trimmed_surface(a.surfaces()[iA], loopsA[iA], proxyB, false, false));

		for (int iB = 0; iB < (int)b.surfaces().size(); ++iB)
			if (touchedB[iB] || insideBInA[iB])
				result.push_back(build_operation_trimmed_surface(b.surfaces()[iB], loopsB[iB], proxyA, true, true));

		return !result.empty();
	}
}

NurbsBoolean::NurbsBoolean()
{
}

bool NurbsBoolean::boolean_union(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const
{
	result.clear();

	const BoundingBox bboxA = compute_bbox(a);
	const BoundingBox bboxB = compute_bbox(b);

	switch (classify_bbox_relation(bboxA, bboxB))
	{
	case BBoxRelation::EmptyA:
		result.append(b);
		return true;
	case BBoxRelation::EmptyB:
		result.append(a);
		return true;
	case BBoxRelation::Disjoint:
		result.append(a);
		result.append(b);
		return true;
	case BBoxRelation::AContainsB:
		result.append(a);
		return true;
	case BBoxRelation::BContainsA:
		result.append(b);
		return true;
	case BBoxRelation::PartialOverlap:
	default:
		return false;
	}
}

bool NurbsBoolean::boolean_intersection(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const
{
	result.clear();

	const BoundingBox bboxA = compute_bbox(a);
	const BoundingBox bboxB = compute_bbox(b);

	switch (classify_bbox_relation(bboxA, bboxB))
	{
	case BBoxRelation::EmptyA:
	case BBoxRelation::EmptyB:
	case BBoxRelation::Disjoint:
		return true;
	case BBoxRelation::AContainsB:
		result.append(b);
		return true;
	case BBoxRelation::BContainsA:
		result.append(a);
		return true;
	case BBoxRelation::PartialOverlap:
	default:
		return false;
	}
}

//returns true if the result is fully determined by the bounding box relation, false if an actual boolean operation is needed
bool NurbsBoolean::boolean_difference_bbox(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result) const
{
	result.clear();

	const BoundingBox bboxA = compute_bbox(a);
	const BoundingBox bboxB = compute_bbox(b);

	switch (classify_bbox_relation(bboxA, bboxB))
	{
	case BBoxRelation::EmptyA:
		return true;
	case BBoxRelation::EmptyB:
	case BBoxRelation::Disjoint:
		result.append(a);
		return true;
	case BBoxRelation::AContainsB:
		// Cannot determine a trimmed difference from bounding boxes when B is contained inside A.
		return false;
	case BBoxRelation::BContainsA:
		return true;
	case BBoxRelation::PartialOverlap:
	default:
		return false;
	}
}

//////////////////////////////////////////////////////////////////////////////////
bool NurbsBoolean::compute_union(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics) const
{
	reset_diagnostics(diagnostics);

	NurbsSolid solidResult;
	if (boolean_union(a, b, solidResult))
	{
		result = solidResult;
		return true;
	}

	std::vector<NurbsTrimmedSurface> trimmedResult;
	if (!build_partial_overlap_exact_trimmed(a, b, ExactTrimmedOperation::Union, trimmedResult, diagnostics))
		return false;
	result.clear();
	for (const auto& ts : trimmedResult)
		result.add_surface(ts);
	return true;
}
//////////////////////////////////////////////////////////////////////////////////
bool NurbsBoolean::compute_intersection(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics) const
{
	reset_diagnostics(diagnostics);

	NurbsSolid solidResult;
	if (boolean_intersection(a, b, solidResult))
	{
		result = solidResult;
		return true;
	}

	std::vector<NurbsTrimmedSurface> trimmedResult;
	if (!build_partial_overlap_exact_trimmed(a, b, ExactTrimmedOperation::Intersection, trimmedResult, diagnostics))
		return false;
	result.clear();
	for (const auto& ts : trimmedResult)
		result.add_surface(ts);
	return true;
}
//////////////////////////////////////////////////////////////////////////////////
bool NurbsBoolean::compute_difference(const NurbsSolid& a, const NurbsSolid& b, NurbsSolid& result, NurbsIntersectionResult* diagnostics) const
{
	reset_diagnostics(diagnostics);

	NurbsSolid solidResult;
	if (boolean_difference_bbox(a, b, solidResult))
	{
		result = solidResult;
		return true;
	}

	std::vector<NurbsTrimmedSurface> trimmedResult;
	if (!build_partial_overlap_exact_trimmed(a, b, ExactTrimmedOperation::Difference, trimmedResult, diagnostics))
		return false;
	result.clear();
	for (const auto& ts : trimmedResult)
		result.add_surface(ts);
	return true;
}

bool NurbsBoolean::compute_difference(const NurbsSolid& a, const NurbsSolid& b, std::vector<NurbsTrimmedSurface>& result, NurbsIntersectionResult* diagnostics) const
{
	reset_diagnostics(diagnostics);

	NurbsSolid solidResult;
	if (boolean_difference_bbox(a, b, solidResult))
	{
		result.clear();
		result.reserve(solidResult.surfaces().size());
		for (const auto& s : solidResult.surfaces())
			result.emplace_back(s);
		return true;
	}

	return build_partial_overlap_exact_trimmed(a, b, ExactTrimmedOperation::Difference, result, diagnostics);
}
//////////////////////////////////////////////////////////////////////////////////