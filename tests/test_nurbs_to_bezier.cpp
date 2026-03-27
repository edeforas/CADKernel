#include "NurbsUtil.h"
#include "NurbsFactory.h"
#include "NurbsSurface.h"
#include "BezierSurface.h"
#include "Mesh.h"

#include <iostream>
#include <vector>

void test_nurbs_to_bezier_conversion()
{
    std::cout << std::endl << "test_nurbs_to_bezier_conversion" << std::endl;

    // Create a simple NURBS surface (a plane)
    NurbsSurface nurbs;
    std::vector<Point3> points = {
        Point3(0, 0, 0), Point3(1, 0, 0), Point3(2, 0, 0),
        Point3(0, 1, 0), Point3(1, 1, 0.5), Point3(2, 1, 0.2),
        Point3(0, 2, 0), Point3(1, 2, 0.1), Point3(2, 2, 0)
    };

    nurbs.set_degree(2, 2);
    nurbs.set_knots_u({0, 0, 0, 1, 1, 1});
    nurbs.set_knots_v({0, 0, 0, 1, 1, 1});
    nurbs.set_points(points, 3, 3);

    // Convert to Bezier patches
    std::vector<BezierSurface> bezierPatches;
    NurbsUtil::to_bezier_patches(nurbs, bezierPatches);

    std::cout << "NURBS surface converted to " << bezierPatches.size() << " Bezier patches" << std::endl;

    // Test that we got patches
    if (bezierPatches.size() > 0) {
        const BezierSurface& firstPatch = bezierPatches[0];
        std::cout << "First patch: degree (" << firstPatch.degree_u() << ", " << firstPatch.degree_v() << ")" << std::endl;
        std::cout << "Control points: " << firstPatch.nb_control_points() << std::endl;

        // Test evaluation
        Point3 p = firstPatch.evaluate(0.5, 0.5);
        std::cout << "Evaluation at (0.5, 0.5): (" << p.x() << ", " << p.y() << ", " << p.z() << ")" << std::endl;

        // Convert to mesh for visualization
        Mesh mesh;
        for (const auto& patch : bezierPatches) {
            // Add patch control points to mesh for visualization
            for (const auto& cp : patch.control_points()) {
                mesh.add_vertex(cp);
            }
        }
        std::cout << "Created mesh with " << mesh.nb_vertices() << " vertices" << std::endl;
    }

    std::cout << "Test completed successfully!" << std::endl;
}

int main()
{
    test_nurbs_to_bezier_conversion();
    return 0;
}