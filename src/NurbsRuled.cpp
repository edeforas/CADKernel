#include "NurbsRuled.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"

#include <algorithm>

///////////////////////////////////////////////////////////////////////////
NurbsRuled::NurbsRuled()
{}
///////////////////////////////////////////////////////////////////////////
NurbsRuled::~NurbsRuled()
{}
///////////////////////////////////////////////////////////////////////////
bool NurbsRuled::create_ruled(const NurbsCurve& nc1, const NurbsCurve& nc2, NurbsSurface& ns) const
{
	ns.clear();

	// Make curves compatible by ensuring they have the same number of control points
	// and knot vectors. Use degree elevation and knot insertion as needed.

	NurbsCurve c1 = nc1;
	NurbsCurve c2 = nc2;

	// Ensure both curves have the same degree
	int maxDegree = std::max(c1.degree(), c2.degree());
	while (c1.degree() < maxDegree) {
		if (!c1.degree_elevation()) return false;
	}
	while (c2.degree() < maxDegree) {
		if (!c2.degree_elevation()) return false;
	}

	// Ensure both curves have the same number of control points
	int maxPoints = std::max(c1.nb_points(), c2.nb_points());
	while (c1.nb_points() < maxPoints) {
		// Insert knots to increase control points
		if (c1.nb_points() < maxPoints) {
			// Insert knot at midpoint of existing spans
			for (size_t i = c1.degree(); i < c1.knots().size() - c1.degree(); ++i) {
				if (c1.knots()[i] < c1.knots()[i+1]) {
					double midKnot = (c1.knots()[i] + c1.knots()[i+1]) * 0.5;
					c1.insert_knot(midKnot);
					if (c1.nb_points() >= maxPoints) break;
				}
			}
		}
	}
	while (c2.nb_points() < maxPoints) {
		for (size_t i = c2.degree(); i < c2.knots().size() - c2.degree(); ++i) {
			if (c2.knots()[i] < c2.knots()[i+1]) {
				double midKnot = (c2.knots()[i] + c2.knots()[i+1]) * 0.5;
				c2.insert_knot(midKnot);
				if (c2.nb_points() >= maxPoints) break;
			}
		}
	}

	// Now both curves should have compatible structure
	if (c1.nb_points() != c2.nb_points()) {
		// Final fallback: use the minimum number of points
		int minPoints = std::min(c1.nb_points(), c2.nb_points());
		// This is a simplification - in practice, we'd need more sophisticated compatibility
		return false; // For now, require exact compatibility
	}

	std::vector<Point3> vp;
	std::vector<double> vw;

	vp = c1.points();
	vp.insert(vp.end(), c2.points().begin(), c2.points().end());

	vw = c1.weights();
	vw.insert(vw.end(), c2.weights().begin(), c2.weights().end());

	std::vector<double> knotsV = { 0., 0., 1., 1. };

	ns.set_degree(c1.degree(), 1);
	ns.set_points(vp, c1.nb_points(), 2);
	ns.set_knots_u(c1.knots());
	ns.set_knots_v(knotsV);
	ns.set_weights(vw);
	ns.set_closed_u(c1.is_closed() && c2.is_closed());
	ns.set_closed_v(false);

	return true;
}
///////////////////////////////////////////////////////////////////////////