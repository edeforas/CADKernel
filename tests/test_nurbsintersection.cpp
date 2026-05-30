#include "NurbsIntersection.h"
#include "NurbsFactory.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

void test(bool b, string s = "")
{
    if (!b)
    {
        cout << "condition: [" << s << "] is not realized!" << endl;
        exit(-1);
    }
}

void test_plane_plane()
{
    cout << "\n test_plane_plane" << endl;
    NurbsSurface s1, s2;

    // large quad in XY plane z=0
    NurbsFactory::create_quad(Point3(-10., -10., 0.), Point3(10., -10., 0.), Point3(10., 10., 0.), Point3(-10., 10., 0.), s1);
    // large quad in XZ plane y=0
    NurbsFactory::create_quad(Point3(-10., 0., -10.), Point3(10., 0., -10.), Point3(10., 0., 10.), Point3(-10., 0., 10.), s2);

    NurbsIntersectionResult res;
    IntersectionOptions opt;
    opt.geomTol = 1e-6;
    opt.seedSamplingRes = 32;

    compute_surface_intersection(s1, s2, res, opt);

    test(res.curves.size() > 0, "no intersection curves found for plane-plane");
    test(res.curves[0].samples.size() > 0, "no samples in found curve");
}

void test_plane_sphere()
{
    cout << "\n test_plane_sphere" << endl;
    NurbsSurface sPlane, sSphere;
    // plane z=0
    NurbsFactory::create_quad(Point3(-6., -6., 0.), Point3(6., -6., 0.), Point3(6., 6., 0.), Point3(-6., 6., 0.), sPlane);
    // sphere radius 5
    NurbsFactory::create_sphere(5.0, sSphere);

    NurbsIntersectionResult res;
    IntersectionOptions opt;
    opt.geomTol = 1e-4; // sphere-surface conversions may be coarser
    opt.seedSamplingRes = 48;

    compute_surface_intersection(sPlane, sSphere, res, opt);

    test(res.curves.size() > 0, "no intersection curves found for plane-sphere");
}

int main()
{
    test_plane_plane();
    test_plane_sphere();
    cout << "nurbs intersection tests passed" << endl;
    return 0;
}
