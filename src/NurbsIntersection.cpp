#include "NurbsIntersection.h"

#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>

namespace {
    bool find_surface_intersection_points(const NurbsSurface& s1, const NurbsSurface& s2,
                                        std::vector<Point3>& intersectionPoints,
                                        std::vector<std::pair<double, double>>& params1,
                                        std::vector<std::pair<double, double>>& params2);

    bool connect_intersection_points(const std::vector<Point3>& points,
                                   const std::vector<std::pair<double, double>>& params1,
                                   const std::vector<std::pair<double, double>>& params2,
                                   NurbsIntersectionCurve& curve);

    double point_surface_distance(const Point3& point, const NurbsSurface& surface,
                                double& u, double& v, double tol = 1e-6);

    bool march_along_curve(const NurbsSurface& s1, const NurbsSurface& s2,
                          double u1_start, double v1_start, double u2_start, double v2_start,
                          std::vector<Point3>& points,
                          std::vector<std::pair<double, double>>& params1,
                          std::vector<std::pair<double, double>>& params2,
                          double step_size = 0.01, int max_steps = 1000);

    bool find_surface_intersection_points(const NurbsSurface& s1, const NurbsSurface& s2,
                                        std::vector<Point3>& intersectionPoints,
                                        std::vector<std::pair<double, double>>& params1,
                                        std::vector<std::pair<double, double>>& params2)
    {
        intersectionPoints.clear();
        params1.clear();
        params2.clear();

        // Simple grid sampling approach for finding intersection points
        // In a full implementation, this would use more sophisticated methods

        const int grid_size = 20; // Sample points per surface
        const double tol = 1e-6;

        for (int i = 0; i <= grid_size; ++i) {
            for (int j = 0; j <= grid_size; ++j) {
                double u1 = static_cast<double>(i) / grid_size;
                double v1 = static_cast<double>(j) / grid_size;

                Point3 p1;
                s1.evaluate_clamped(u1, v1, p1);

                // Project p1 onto s2
                double u2_proj, v2_proj;
                double dist = point_surface_distance(p1, s2, u2_proj, v2_proj, tol);

                if (dist < tol) {
                    // Check if this point is already found (within tolerance)
                    bool duplicate = false;
                    for (size_t k = 0; k < intersectionPoints.size(); ++k) {
                        if (intersectionPoints[k].distance_square(p1) < tol * tol) {
                            duplicate = true;
                            break;
                        }
                    }

                    if (!duplicate) {
                        intersectionPoints.push_back(p1);
                        params1.emplace_back(u1, v1);
                        params2.emplace_back(u2_proj, v2_proj);
                    }
                }
            }
        }

        return !intersectionPoints.empty();
    }

    bool connect_intersection_points(const std::vector<Point3>& points,
                                   const std::vector<std::pair<double, double>>& params1,
                                   const std::vector<std::pair<double, double>>& params2,
                                   NurbsIntersectionCurve& curve)
    {
        if (points.empty()) return false;

        curve.samples.clear();
        curve.closed = false;

        // For now, just add all points as samples
        // In a full implementation, this would sort and connect them properly
        for (size_t i = 0; i < points.size(); ++i) {
            NurbsIntersectionSample sample;
            sample.point = points[i];
            sample.uA = params1[i].first;
            sample.vA = params1[i].second;
            sample.uB = params2[i].first;
            sample.vB = params2[i].second;
            curve.samples.push_back(sample);
        }

        return true;
    }

    double point_surface_distance(const Point3& point, const NurbsSurface& surface,
                                double& u, double& v, double tol)
    {
        // Simple projection using surface evaluation
        // In a full implementation, this would use Newton iteration or similar

        double min_dist = std::numeric_limits<double>::max();
        double best_u = 0.5, best_v = 0.5;

        const int samples = 10;
        for (int i = 0; i <= samples; ++i) {
            for (int j = 0; j <= samples; ++j) {
                double test_u = static_cast<double>(i) / samples;
                double test_v = static_cast<double>(j) / samples;

                Point3 p;
                surface.evaluate_clamped(test_u, test_v, p);

                double dist = point.distance_square(p);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_u = test_u;
                    best_v = test_v;
                }
            }
        }

        u = best_u;
        v = best_v;
        return std::sqrt(min_dist);
    }

    bool march_along_curve(const NurbsSurface& s1, const NurbsSurface& s2,
                          double u1_start, double v1_start, double u2_start, double v2_start,
                          std::vector<Point3>& points,
                          std::vector<std::pair<double, double>>& params1,
                          std::vector<std::pair<double, double>>& params2,
                          double step_size, int max_steps)
    {
        // Basic marching algorithm - in a full implementation this would be more sophisticated
        double u1 = u1_start, v1 = v1_start;
        double u2 = u2_start, v2 = v2_start;

        points.push_back(Point3()); // Placeholder
        params1.emplace_back(u1, v1);
        params2.emplace_back(u2, v2);

        // This is a simplified implementation
        return true;
    }
}

// Main intersection function implementation
void compute_surface_intersection(const NurbsSurface& surface1, const NurbsSurface& surface2,
                                NurbsIntersectionResult& result)
{
    result.curves.clear();
    result.hasPartialOverlap = false;

    std::vector<Point3> intersectionPoints;
    std::vector<std::pair<double, double>> params1, params2;

    if (find_surface_intersection_points(surface1, surface2, intersectionPoints, params1, params2)) {
        NurbsIntersectionCurve curve;
        if (connect_intersection_points(intersectionPoints, params1, params2, curve)) {
            result.curves.push_back(curve);
        }
    }

    // For now, assume no partial overlap
    result.hasPartialOverlap = false;
}