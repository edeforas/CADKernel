#ifndef _MeshFactory_
#define _MeshFactory_

#include "Mesh.h"
#include <vector>

namespace MeshFactory
{
    void create_tetrahedron(double dSize, Mesh& m);
    void create_box(double xSize, double ySize, double zSize, Mesh& m);
    void create_octahedron(double dSize, Mesh& m);
    void create_dodecahedron(double dSize, Mesh& m);
    void create_icosahedron(double dSize, Mesh& m);
    void create_cylinder(double dHeight, double dRadius,int iNbSegments, Mesh& m);
    void create_sphere_geodesic(double dRadius, int iNbSegments, Mesh& m);
    void create_sphere_uv(double dRadius, int iNbSegments, Mesh& m);
    void create_torus(double dMajorRadius, double dMinorRadius, int iNbSegments, Mesh& m);
    void create_revolve(const std::vector<Point3>& profile, int iNbSegments, Mesh& m);
    void create_revolve(const std::vector<Point3>& profile, int iNbSegments, bool bCapEnds, Mesh& m);
}
#endif