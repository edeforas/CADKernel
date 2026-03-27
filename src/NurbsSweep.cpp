#include "NurbsSweep.h"

#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <cmath>
#include <vector>


bool NurbsSweep::sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface, bool perpendicular)
{
    surface.clear();

    int nU = profile.nb_points();
    int nV = path.nb_points();

    if (nU <= 0 || nV <= 0)
        return false;

    const std::vector<Point3>& profilePoints = profile.points();
    const std::vector<Point3>& pathPoints = path.points();
    if ((int)profilePoints.size() != nU || (int)pathPoints.size() != nV)
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
    std::vector<double> pathWeights = NurbsUtil::build_safe_weights(path.weights(), nV);

    surface.set_degree(degreeU, degreeV);
    surface.set_knots_u(knotsU);
    surface.set_knots_v(knotsV);

    std::vector<Point3> surfacePoints;
    std::vector<double> surfaceWeights;
    surfacePoints.reserve(nU * nV);
    surfaceWeights.reserve(nU * nV);

    Point3 pathStart = pathPoints[0];
    pathStart.sanitize();
    for (int j = 0; j < nV; ++j)
    {
        Point3 currentPathPoint = pathPoints[j];
        currentPathPoint.sanitize();
        double pathWeight = pathWeights[j];
        Point3 offset = currentPathPoint - pathStart;

        // If perpendicular, compute rotation
        Point3 rotationAxis(0, 0, 0);
        double rotationAngle = 0.0;
        if (perpendicular) {
            // Compute tangent using finite differences
            Point3 tangent(0, 0, 0);
            if (j == 0 && nV > 1) {
                tangent = pathPoints[1] - pathPoints[0];
            } else if (j == nV - 1 && nV > 1) {
                tangent = pathPoints[nV - 1] - pathPoints[nV - 2];
            } else if (nV > 2) {
                tangent = (pathPoints[j + 1] - pathPoints[j - 1]) * 0.5;
            }
            if (tangent.norm_square() > 1e-12) {
                tangent.normalize();
                // Compute new normal: perpendicular to tangent
                Point3 zAxis(0, 0, 1);
                Point3 newNormal = zAxis.cross_product(tangent);
                double newNormalNormSq = newNormal.norm_square();
                if (newNormalNormSq > 1e-12) {
                    newNormal.normalize();
                    // Rotation axis: zAxis × newNormal
                    rotationAxis = zAxis.cross_product(newNormal);
                    double axisNormSq = rotationAxis.norm_square();
                    if (axisNormSq > 1e-12) {
                        rotationAxis.normalize();
                        // Angle: acos(zAxis • newNormal)
                        double cosAngle = zAxis.dot_product(newNormal);
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        rotationAngle = acos(cosAngle);
                    }
                }
                // If newNormal is zero (tangent parallel to z-axis), no rotation needed
            }
        }

        for (int i = 0; i < nU; ++i)
        {
            Point3 profilePoint = profilePoints[i];
            profilePoint.sanitize();
            double profileWeight = profileWeights[i];

            // Apply rotation if perpendicular
            if (perpendicular && rotationAxis.norm_square() > 1e-12 && fabs(rotationAngle) > 1e-12) {
                // Rodrigues' rotation formula
                double cosA = cos(rotationAngle);
                double sinA = sin(rotationAngle);
                Point3 axis = rotationAxis; // Already normalized
                Point3 rotated = profilePoint * cosA +
                                axis.cross_product(profilePoint) * sinA +
                                axis * (axis.dot_product(profilePoint) * (1 - cosA));
                profilePoint = rotated;
                profilePoint.sanitize(); // Ensure no NaN/inf values
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