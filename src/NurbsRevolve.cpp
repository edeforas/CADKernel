#include "NurbsRevolve.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "NurbsBasis.h"

#include <cmath>
#include <algorithm>
#include <vector>

namespace
{
	constexpr double kPi = NurbsBasis::Pi;

	Point3 rotate_about_z(const Point3& p, double a)
	{
		const double c = std::cos(a);
		const double s = std::sin(a);
		return Point3(c * p.x() - s * p.y(), s * p.x() + c * p.y(), p.z());
	}
}

///////////////////////////////////////////////////////////////////////////
NurbsRevolve::NurbsRevolve()
{
}
///////////////////////////////////////////////////////////////////////////
NurbsRevolve::~NurbsRevolve()
{
}
///////////////////////////////////////////////////////////////////////////
bool NurbsRevolve::revolve(const NurbsCurve& nc, NurbsSurface& ns) const
{
	return revolve(nc, 2. * kPi, ns);
}
///////////////////////////////////////////////////////////////////////////
bool NurbsRevolve::revolve(const NurbsCurve& nc, double dAngleRad, NurbsSurface& ns) const
{
	ns.clear();

	const auto& refPoints = nc.points();
	if (refPoints.empty())
		return false;

	const auto& refWeights = nc.weights();
	if ((int)refWeights.size() != (int)refPoints.size())
		return false;

	if (!std::isfinite(dAngleRad))
		return false;

	const double absAngle = std::fabs(dAngleRad);
	if (absAngle < 1.e-12)
		return false;

	double sweep = dAngleRad;
	bool bClosed = false;
	if (absAngle >= (2. * kPi - 1.e-9))
	{
		sweep = (dAngleRad >= 0.) ? 2. * kPi : -2. * kPi;
		bClosed = true;
	}

	const int nbSegments = std::max(1, (int)std::ceil(std::fabs(sweep) / (kPi / 2.)));
	const int nbCtrlU = 2 * nbSegments + 1;
	const double dSeg = sweep / (double)nbSegments;
	const double localMidWeight = std::cos(std::fabs(dSeg) * 0.5);

	std::vector<Point3> pc;
	std::vector<double> wc;
	pc.reserve(nbCtrlU * refPoints.size());
	wc.reserve(nbCtrlU * refPoints.size());

	for (int i = 0; i < (int)refPoints.size(); ++i)
	{
		const Point3 pBase = refPoints[i];
		const double wBase = refWeights[i];

		for (int iSeg = 0; iSeg < nbSegments; ++iSeg)
		{
			const double a0 = iSeg * dSeg;
			const double a1 = (iSeg + 1) * dSeg;

			const Point3 p0 = rotate_about_z(pBase, a0);
			const Point3 p2 = rotate_about_z(pBase, a1);

			if (iSeg == 0)
			{
				pc.push_back(p0);
				wc.push_back(wBase);
			}

			const double r2 = p0.x() * p0.x() + p0.y() * p0.y();
			if (r2 <= 1.e-20)
			{
				pc.push_back(p0);
				wc.push_back(wBase * localMidWeight);
				pc.push_back(p2);
				wc.push_back(wBase);
				continue;
			}

			const Point3 t0(-p0.y(), p0.x(), 0.);
			const Point3 t2(-p2.y(), p2.x(), 0.);

			double den = t0.x() * t2.y() - t0.y() * t2.x();
			if (std::fabs(den) < 1.e-15)
			{
				const Point3 pm = rotate_about_z(pBase, 0.5 * (a0 + a1));
				pc.push_back(pm);
				wc.push_back(wBase * localMidWeight);
				pc.push_back(p2);
				wc.push_back(wBase);
				continue;
			}

			const Point3 d = p2 - p0;
			const double alpha = (d.x() * t2.y() - d.y() * t2.x()) / den;
			Point3 p1 = p0 + t0 * alpha;
			p1.z() = p0.z();

			pc.push_back(p1);
			wc.push_back(wBase * localMidWeight);
			pc.push_back(p2);
			wc.push_back(wBase);
		}
	}

	ns.set_degree(2, nc.degree());
	ns.set_points(pc, nbCtrlU, nc.nb_points());
	ns.set_knots_u(NurbsUtil::build_segmented_quadratic_knots(nbSegments));
	ns.set_knots_v(nc.knots());
	ns.set_weights(wc);

	ns.set_closed_u(bClosed);
	ns.set_closed_v(nc.is_closed());

	return true;
}
///////////////////////////////////////////////////////////////////////////
bool NurbsRevolve::revolve_deg(const NurbsCurve& nc, double dAngleDeg, NurbsSurface& ns) const
{
	if (!std::isfinite(dAngleDeg))
		return false;

	return revolve(nc, dAngleDeg * kPi / 180., ns);
}
///////////////////////////////////////////////////////////////////////////