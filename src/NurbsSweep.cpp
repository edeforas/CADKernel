#include "NurbsSweep.h"

#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <cmath>
#include <vector>


bool NurbsSweep::sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface)
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

        for (int i = 0; i < nU; ++i)
        {
            Point3 profilePoint = profilePoints[i];
            profilePoint.sanitize();
            double profileWeight = profileWeights[i];

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