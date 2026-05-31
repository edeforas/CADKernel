#include "NurbsSweep.h"

#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsUtil.h"
#include "BoundingBox.h"
#include "NurbsFactory.h"
#include "NurbsSolid.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
	// Helper function to compute an orthonormal frame aligned with a given normal direction
	void compute_local_frame(const Point3& normal, Point3& localX, Point3& localY, Point3& localZ)
	{
		localZ = normal;
		if (localZ.norm_square() > 1e-12) {
			localZ.normalize();
		} else {
			localZ = Point3(0, 0, 1);
		}
		
		// Create orthonormal basis with localZ as the normal
		Point3 zAxis(0, 0, 1);
		localX = Point3(1, 0, 0);
		localY = Point3(0, 1, 0);
		
		if (fabs(localZ.dot_product(zAxis)) < 0.99) {
			localX = zAxis.cross_product(localZ);
			if (localX.norm_square() < 1e-12) {
				localX = Point3(1, 0, 0);
			}
			localX.normalize();
			localY = localZ.cross_product(localX);
			localY.normalize();
		} else {
			localX = Point3(1, 0, 0);
			if (fabs(localX.dot_product(localZ)) > 0.99) {
				localX = Point3(0, 1, 0);
			}
			localX = localX - localZ * localX.dot_product(localZ);
			if (localX.norm_square() > 1e-12) {
				localX.normalize();
			} else {
				localX = Point3(1, 0, 0);
			}
			localY = localZ.cross_product(localX);
			localY.normalize();
		}
	}

	// Helper function to create a disk cap from the profile curve at a given position
	void create_disk_cap_from_profile(const NurbsCurve& profile, const Point3& position, const Point3& normal, NurbsSurface& cap)
	{
		cap.clear();
		
		// Create a filled surface from the profile curve
		NurbsUtil::create_surface_from_curve(profile, cap);
		
		// Compute local frame aligned with the normal direction
		Point3 localX, localY, localZ;
		compute_local_frame(normal, localX, localY, localZ);
		
		// Transform all cap points: rotate then translate
		for (auto& p : cap.points()) {
			Point3 rotated(
				p.x() * localX.x() + p.y() * localY.x(),
				p.x() * localX.y() + p.y() * localY.y(),
				p.x() * localX.z() + p.y() * localY.z()
			);
			p = rotated + position;
		}
	}
}


bool NurbsSweep::sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface, bool perpendicular)
{
    surface.clear();

    int nU = profile.nb_points();
    int nV = path.nb_points();

    if (nU <= 0 || nV <= 0)
        return false;

    const std::vector<Point3>& profilePoints = profile.points();
    if ((int)profilePoints.size() != nU)
        return false;

    int degreeU = profile.degree();
    int degreeV = path.degree();
    degreeU = std::max(0, std::min(degreeU, nU - 1));
    degreeV = std::max(0, std::min(degreeV, nV - 1));

    std::vector<double> knotsU = profile.knots();
    if (!NurbsCurveUtil::has_valid_knot_vector(knotsU, degreeU, nU))
        knotsU = NurbsCurveUtil::build_clamped_uniform_knots(nU, degreeU);

    std::vector<double> knotsV = path.knots();
    if (!NurbsCurveUtil::has_valid_knot_vector(knotsV, degreeV, nV))
        knotsV = NurbsCurveUtil::build_clamped_uniform_knots(nV, degreeV);

    std::vector<double> profileWeights = NurbsUtil::build_safe_weights(profile.weights(), nU);
    std::vector<double> pathControlWeights = NurbsUtil::build_safe_weights(path.weights(), nV);
    std::vector<double> pathWeights(nV, 1.0);

    NurbsCurve safePath = path;
    safePath.set_degree(degreeV);
    safePath.set_points(path.points());
    safePath.set_knots(knotsV);
    safePath.set_weights(pathControlWeights);

    std::vector<Point3> pathPoints(nV);
    for (int j = 0; j < nV; ++j)
    {
        double t = (nV > 1) ? (double)j / (double)(nV - 1) : 0.0;
        safePath.evaluate(t, pathPoints[j]);
    }

    surface.set_degree(degreeU, degreeV);
    surface.set_knots_u(knotsU);
    surface.set_knots_v(knotsV);

    std::vector<Point3> surfacePoints;
    std::vector<double> surfaceWeights;
    surfacePoints.reserve(nU * nV);
    surfaceWeights.reserve(nU * nV);

    Point3 pathStart = pathPoints[0];
    pathStart.sanitize();

    Point3 prevTangent(0, 0, 0);
    Point3 prevNormal(0, 0, 0);
    Point3 prevBinormal(0, 0, 0);
    bool hasPrevFrame = false;

    for (int j = 0; j < nV; ++j)
    {
        Point3 currentPathPoint = pathPoints[j];
        currentPathPoint.sanitize();
        double pathWeight = pathWeights[j];
        Point3 offset = currentPathPoint - pathStart;

        Point3 tangent(0, 0, 0);
        if (j == 0 && nV > 1) {
            tangent = pathPoints[1] - pathPoints[0];
        } else if (j == nV - 1 && nV > 1) {
            tangent = pathPoints[nV - 1] - pathPoints[nV - 2];
        } else if (nV > 2) {
            tangent = (pathPoints[j + 1] - pathPoints[j - 1]) * 0.5;
        }

        Point3 localX(1, 0, 0);
        Point3 localY(0, 1, 0);
        Point3 localZ(0, 0, 1);
        if (perpendicular && tangent.norm_square() > 1e-12) {
            tangent.normalize();

            Point3 normal(0, 0, 0);
            Point3 binormal(0, 0, 0);
            Point3 zAxis(0, 0, 1);
            Point3 xAxis(1, 0, 0);
            Point3 yAxis(0, 1, 0);

            if (!hasPrevFrame) {
                if (fabs(tangent.dot_product(zAxis)) < 0.99) {
                    normal = zAxis - tangent * tangent.dot_product(zAxis);
                } else {
                    normal = xAxis - tangent * tangent.dot_product(xAxis);
                }
                if (normal.norm_square() < 1e-12) {
                    normal = xAxis;
                }
                normal.normalize();
                binormal = tangent.cross_product(normal);
                if (binormal.norm_square() < 1e-12) {
                    binormal = yAxis - tangent * tangent.dot_product(yAxis);
                }
                binormal.normalize();
                hasPrevFrame = true;
            } else {
                normal = prevNormal - tangent * tangent.dot_product(prevNormal);
                if (normal.norm_square() < 1e-12) {
                    normal = prevBinormal.cross_product(tangent);
                }
                if (normal.norm_square() < 1e-12) {
                    normal = xAxis - tangent * tangent.dot_product(xAxis);
                }
                normal.normalize();
                binormal = tangent.cross_product(normal);
                if (binormal.norm_square() < 1e-12) {
                    binormal = yAxis - tangent * tangent.dot_product(yAxis);
                }
                binormal.normalize();
            }

            prevTangent = tangent;
            prevNormal = normal;
            prevBinormal = binormal;

            localX = normal;
            localY = binormal;
            localZ = tangent;
        }

        for (int i = 0; i < nU; ++i)
        {
            Point3 profilePoint = profilePoints[i];
            profilePoint.sanitize();
            double profileWeight = profileWeights[i];

            if (perpendicular && tangent.norm_square() > 1e-12) {
                Point3 rotated(
                    profilePoint.x() * localX.x() + profilePoint.y() * localY.x() + profilePoint.z() * localZ.x(),
                    profilePoint.x() * localX.y() + profilePoint.y() * localY.y() + profilePoint.z() * localZ.y(),
                    profilePoint.x() * localX.z() + profilePoint.y() * localY.z() + profilePoint.z() * localZ.z());
                profilePoint = rotated;
                profilePoint.sanitize();
            }

            surfacePoints.push_back(profilePoint + offset);
            surfaceWeights.push_back(profileWeight * pathWeight);
        }
    }

    surface.set_points(surfacePoints, nU, nV);
    surface.set_weights(surfaceWeights);
    surface.set_closed_u(profile.is_closed());
    surface.set_closed_v(path.is_closed());

    return true;
}

bool NurbsSweep::sweep_solid(const NurbsCurve& profile, const NurbsCurve& path, NurbsSolid& solid, bool perpendicular)
{
	solid.clear();
	
	// Create and add main sweep surface
	NurbsSurface mainSurface;
	if (!sweep(profile, path, mainSurface, perpendicular)) {
		return false;
	}
	solid.add_surface(mainSurface);
	
	// Get path endpoints and compute tangent vectors
	const std::vector<Point3>& pathPoints = path.points();
	Point3 pathStart = pathPoints[0];
	Point3 pathEnd = pathPoints[pathPoints.size() - 1];
	
	Point3 startTangent = (pathPoints.size() > 1) ? 
		(pathPoints[1] - pathPoints[0]) : 
		Point3(0, 0, 1);
	startTangent.normalize();
	
	Point3 endTangent = (pathPoints.size() > 1) ? 
		(pathPoints[pathPoints.size() - 1] - pathPoints[pathPoints.size() - 2]) : 
		Point3(0, 0, 1);
	endTangent.normalize();
	
	// Compute profile radius for caps
	const std::vector<Point3>& profilePoints = profile.points();
	const std::vector<double>& profileWeights = profile.weights();
	
	Point3 profileCenter(0, 0, 0);
	double totalWeight = 0.0;
	for (int i = 0; i < (int)profilePoints.size(); ++i) {
		double w = (profileWeights.size() == profilePoints.size()) ? profileWeights[i] : 1.0;
		profileCenter += profilePoints[i] * w;
		totalWeight += w;
	}
	if (totalWeight > 1e-12) {
		profileCenter /= totalWeight;
	} else if (profilePoints.size() > 0) {
		profileCenter = profilePoints[0];
	}
	
	double profileRadius = 0.0;
	for (const auto& p : profilePoints) {
		profileRadius = std::max(profileRadius, (p - profileCenter).norm());
	}
	profileRadius = std::max(profileRadius, 1e-12);
	
	// Create and add start cap (disk)
	{
		NurbsSurface capSurface;
		Point3 startNormal = startTangent * -1.0;  // Point backwards from start
		create_disk_cap_from_profile(profile, pathStart, startNormal, capSurface);
		solid.add_surface(capSurface);
	}
	
	// Create and add end cap (disk)
	{
		NurbsSurface capSurface;
		create_disk_cap_from_profile(profile, pathEnd, endTangent, capSurface);
		solid.add_surface(capSurface);
	}
	
	return true;
}
