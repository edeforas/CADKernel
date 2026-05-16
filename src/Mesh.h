#ifndef _Mesh_
#define _Mesh_

#include "Geometry.h"
class MeshKernelLinkedTriangles;
class Transform;

///////////////////////////////////////////////////////////////////////////
class Mesh
{
public:
    Mesh();
    Mesh(const Mesh& m);
    virtual ~Mesh();

    void set_color(int iColor);
    int get_color() const;

    Mesh& operator=(const Mesh& m);
    void add_mesh(const Mesh& f);

    int nb_vertices() const;
    int add_vertex(const Point3& vertex);
    int add_vertex(double a, double b, double c);
    void set_vertex(int iVertex, const Point3& vertex);
    void get_vertex(int iVertex, Point3& vertex) const;

    int add_triangle(int iVertex1, int iVertex2, int iVertex3); // return the triangle ID
    int add_triangle(const Point3& p1, const  Point3& p2, const Point3& p3); //helper
    void get_triangle(int iTriangle, Triangle3& f) const;
    void get_triangle_vertices(int iTriangle, int& iVertex1, int& iVertex2, int& iVertex3) const;
    void get_triangle_vertices(int iTriangle, Point3& p1, Point3& p2, Point3& p3) const;
    int nb_triangles() const;
    void unlink_triangle(int iTriangle);
    bool is_triangle_unlinked(int iTriangle) const;
    void get_near_triangles(int iTriangle, int& iT1, int& iT2, int& iT3) const;

    void split_triangle_with_vertex(int iTriangle, const Point3& p);
    void split_triangle_with_vertex(int iTriangle, int iVertex); //create 3 triangles with edges and new vertex, remove iTriangle
    void split_edge_with_vertex(int iTriangle1, int iTriangle2,int iVertex1, int iVertex2, int iVertexSplit); // 4 new triangles added at the end, unlink the two previous faces
    void split_triangle(int iTriangle, const Triangle3& tSplitter);

    void flip_triangle(int iTriangle); //flip triangle so normal is reversed

    void add_quad(int iVertex1, int iVertex2, int iVertex3, int iVertex4);
    void add_quad(const Point3& p1, const  Point3& p2, const Point3& p3, const Point3& p4,bool bOptimSurface = true); //helper

    void add_pentagon(int iVertex1, int iVertex2, int iVertex3, int iVertex4, int iVertex5);      

    void apply_transform(const Transform& t);

    void clear();
    bool empty() const;

private:
	MeshKernelLinkedTriangles* _pKernel;
    int _iColor; //-1 for no custom color
};
///////////////////////////////////////////////////////////////////////////

#endif