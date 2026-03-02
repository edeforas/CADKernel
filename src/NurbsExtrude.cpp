#include "NurbsExtrude.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"

///////////////////////////////////////////////////////////////////////////
NurbsExtrude::NurbsExtrude()
{}
///////////////////////////////////////////////////////////////////////////
NurbsExtrude::~NurbsExtrude()
{}
///////////////////////////////////////////////////////////////////////////
bool NurbsExtrude::extrude(const NurbsCurve& nc, const Point3& direction, NurbsSurface& ns) const
{
	ns.clear();

	//extrude along V

	std::vector<double> weights = nc.weights();
	std::vector<Point3> points= nc.points();
	int iOldSize = points.size();

	//duplicate weights
	for (int i = 0; i < iOldSize; i++)
		weights.push_back(weights[i]);

	//duplicate points + directions
	for (int i = 0; i < iOldSize; i++)
		points.push_back(points[i]+direction);

	// extrude
	int degreeV = 1;
	std::vector<double> knotsV = { 0., 0., 1., 1. };

	ns.set_degree(nc.degree(), degreeV);
	ns.set_knots_u(nc.knots());
	ns.set_knots_v(knotsV);

	ns.set_weights(weights);
	ns.set_points(points, iOldSize,2);

	ns.set_closed_u(nc.is_closed());

	return true;
}
///////////////////////////////////////////////////////////////////////////