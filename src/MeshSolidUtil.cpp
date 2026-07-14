#include "MeshSolidUtil.h"
#include "Mesh.h"

void MeshSolidUtil::generate_faces(const Mesh& m,  MeshSolid& ms, double dToleranceAngleDeg)
{
    vector<Mesh> vSurfaces;

    //group nearby triangles into surfaces
    ms.clear();

    for (int iTriangle = 0; iTriangle < m.nb_triangles(); iTriangle++)
    {
        if (m.is_triangle_unlinked(iTriangle))
            continue;

        Point3 p1, p2, p3, normal;
        m.get_triangle_vertices(iTriangle, p1, p2, p3);
        m.get_normal(iTriangle, normal);

        bool bFoundSurface = false;
        for(int iS=0; iS<vSurfaces.size() &&(!bFoundSurface); iS++)
        {
            Mesh& surface = vSurfaces[iS];

            for (int jTriangle = 0; jTriangle < surface.nb_triangles() && (!bFoundSurface); jTriangle++)
            {
                Point3 normal2;
                surface.get_normal(jTriangle, normal2);                
                double dAngle = normal.angle_with(normal2);
                if ((dAngle < dToleranceAngleDeg))
                    if(surface.common_edge_with(jTriangle, p1, p2, p3))
                    {    
                        surface.add_triangle(p1, p2, p3);
                        bFoundSurface = true;
                        continue;
                    }
            }
        }

        if(!bFoundSurface)
        {
            Mesh newSurface;
            newSurface.add_triangle(p1, p2, p3);
            vSurfaces.push_back(newSurface);
        }
    }

    for (int iS = 0; iS < vSurfaces.size(); iS++)
    {
        Mesh& surface = vSurfaces[iS];
        ms.add_surface(surface);
    }   
}