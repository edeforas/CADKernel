#include "NurbsFitting.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "NurbsRuled.h"
#include "StepFile.h"
#include "OBJFile.h"

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Parametric equations for a Mobius strip
// from https://en.wikipedia.org/wiki/M%C3%B6bius_strip
Point3 mobius_strip_point(double u, double v, double a)
{
    u *= 2.0 * M_PI;
    v *= 2.0 * M_PI;

    double x = (a + cos(u/2.0)*v/2.0)*cos(u);
    double y = (a + cos(u/2.0)*v/2.0)*sin(u);
    double z = sin(u/2.0)*v/2.0;

    return Point3(x, y, z);
}

int main()
{
    cout << "Creating Mobius strip sample" << endl;

    NurbsCurve n1;
    NurbsCurve n2;
    const double a = 6.0;
    const int numSamplesU = 32;  // Number of samples in u direction
    
    vector<Point3> samplesN1, samplesN2;
    for (int i = 0; i < numSamplesU; ++i)
    {
        double uf=i / (numSamplesU - 1.);
        samplesN1.push_back(mobius_strip_point(uf, -1, a));
        samplesN2.push_back(mobius_strip_point(uf, 1, a));
    }
        
    NurbsFitting::fit_curve_least_squares(samplesN1, 3, 8, n1);
    NurbsFitting::fit_curve_least_squares(samplesN2, 3, 8, n2);

    NurbsSurface ns;
    NurbsRuled nr;
    nr.create_ruled(n1, n2, ns);

    cout << "Successfully created Mobius strip NURBS surface" << endl;
    cout << "Surface degree: (" << ns.degree_u() << ", " << ns.degree_v() << ")" << endl;
    cout << "Control points: " << ns.nb_points_u() << " x " << ns.nb_points_v() << endl;

    // Write to STEP file
    StepWriter sw;
    sw.open("sample_mobius_strip.step");
    sw.write(ns);
    sw.close();

    // Convert to mesh and write to OBJ file
    Mesh mesh;
    NurbsUtil::to_mesh(ns, mesh, 16);  

    OBJWriter ow;
    ow.open("sample_mobius_strip.obj");
    ow.write(mesh);
    ow.close();

    cout << "Mobius strip written to sample_mobius_strip.step and sample_mobius_strip.obj" << endl;
    cout << "Mesh has " << mesh.nb_vertices() << " vertices and " << mesh.nb_triangles() << " triangles" << endl;

    return 0;
}