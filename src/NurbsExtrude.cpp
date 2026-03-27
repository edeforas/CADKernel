#include "NurbsExtrude.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsUtil.h"
#include "Transform.h"

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
	std::vector<Point3> points = nc.points();
	int iOldSize = points.size();

	// Create translated points first, then original points
	std::vector<Point3> extrudedPoints;
	std::vector<double> extrudedWeights;
	for (int i = 0; i < iOldSize; i++) {
		extrudedPoints.push_back(points[i] + direction);
		extrudedWeights.push_back(weights[i]);
	}
	for (int i = 0; i < iOldSize; i++) {
		extrudedPoints.push_back(points[i]);
		extrudedWeights.push_back(weights[i]);
	}

	// extrude
	int degreeV = 1;
	std::vector<double> knotsV = { 0., 0., 1., 1. };

	ns.set_degree(nc.degree(), degreeV);
	ns.set_knots_u(nc.knots());
	ns.set_knots_v(knotsV);

	ns.set_weights(extrudedWeights);
	ns.set_points(extrudedPoints, iOldSize, 2);

	ns.set_closed_u(nc.is_closed());

	return true;
}
///////////////////////////////////////////////////////////////////////////
bool NurbsExtrude::extrude_face(NurbsSolid& solid, int faceIndex, const Point3& direction) const
{
	if (faceIndex < 0 || faceIndex >= (int)solid.surfaces().size()) {
		return false;
	}

	const NurbsSurface& face = solid.surfaces()[faceIndex];

	// Create boundary curve by evaluating face edges (at original position)
	std::vector<Point3> boundaryPoints;
	int numSamples = 10; // Samples per edge

	// Bottom edge: v=0, u from 0 to 1
	for (int j = 0; j < numSamples; ++j) {
		double u = (double)j / (numSamples - 1);
		Point3 p;
		face.evaluate(u, 0.0, p);
		boundaryPoints.push_back(p);
	}
	// Right edge: u=1, v from 0 to 1
	for (int j = 1; j < numSamples; ++j) { // Skip first point to avoid duplicate
		double v = (double)j / (numSamples - 1);
		Point3 p;
		face.evaluate(1.0, v, p);
		boundaryPoints.push_back(p);
	}
	// Top edge: v=1, u from 1 to 0
	for (int j = numSamples - 2; j >= 0; --j) { // Skip last, reverse order
		double u = (double)j / (numSamples - 1);
		Point3 p;
		face.evaluate(u, 1.0, p);
		boundaryPoints.push_back(p);
	}
	// Left edge: u=0, v from 1 to 0
	for (int j = numSamples - 2; j >= 0; --j) { // Skip first, include last to close the loop
		double v = (double)j / (numSamples - 1);
		Point3 p;
		face.evaluate(0.0, v, p);
		boundaryPoints.push_back(p);
	}

	// Create curve from boundary points
	NurbsCurve faceCurve;
	NurbsUtil::create_curve_from_points(boundaryPoints, 1, faceCurve); // Use degree 1 for sharp edges

	// Translate the face
	Translation translation(direction);
	solid.surfaces()[faceIndex].apply_transform(translation);

	// Extrude the curve to connect original boundary to translated boundary
	NurbsSurface extrudedSurface;
	bool success = extrude(faceCurve, direction, extrudedSurface);
	if (success) {
		solid.add_surface(extrudedSurface);
	}

	return success;
}
///////////////////////////////////////////////////////////////////////////