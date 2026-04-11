#include "Mesh.h"

#include <cassert>

#include "MeshKernelLinkedTriangles.h"

#include "Transform.h"

#define swap(a,b) { auto tmp=a; a=b; b=tmp; }  
#define rotate(a,b,c) { auto tmp=a; a=b; b=c; c=a; }  

///////////////////////////////////////////////////////////////////////////
Mesh::Mesh()
{ 
	_pKernel = new MeshKernelLinkedTriangles;
	_iColor = -1;
}

Mesh::Mesh(const Mesh& m)
{
	this->operator=(m);
}

Mesh::~Mesh()
{ 
	delete _pKernel;
}

void Mesh::set_color(int iColor)
{
	_iColor = iColor;
}

int Mesh::get_color() const
{
	return _iColor;
}

int Mesh::nb_triangles() const
{
	return _pKernel->nb_triangles();
}

bool Mesh::is_triangle_unlinked(int iTriangle) const
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	return _pKernel->is_triangle_unlinked(iTriangle);
}
void Mesh::get_near_triangles(int iTriangle, int& iT1, int& iT2, int& iT3) const
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	return _pKernel->get_near_triangles(iTriangle,iT1,iT2,iT3);
}
///////////////////////////////////////////////////////////////////////////
void Mesh::get_triangle(int iTriangle, Triangle3& f) const
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	int iP1, iP2, iP3;
	get_triangle_vertices(iTriangle, iP1, iP2, iP3);
	
	Point3 v1,v2,v3;
	_pKernel->get_vertex(iP1,v1);
	_pKernel->get_vertex(iP2,v2);
	_pKernel->get_vertex(iP3,v3);
	
	f.set_p1(v1);
	f.set_p2(v2);
	f.set_p3(v3);
}

void Mesh::get_triangle_vertices(int iTriangle, int& iVertex1, int& iVertex2, int& iVertex3) const
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	_pKernel->get_triangle_vertices(iTriangle,iVertex1,iVertex2,iVertex3);
}

void Mesh::get_triangle_vertices(int iTriangle, Point3& p1, Point3& p2, Point3& p3) const
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	int iP1, iP2, iP3;
	get_triangle_vertices(iTriangle, iP1, iP2, iP3);

	_pKernel->get_vertex(iP1, p1);
	_pKernel->get_vertex(iP2, p2);
	_pKernel->get_vertex(iP3, p3);
}

void Mesh::unlink_triangle(int iTriangle)
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	_pKernel->unlink_triangle(iTriangle);
}

Mesh& Mesh::operator=(const Mesh& m)
{
	if (this == &m)
		return *this;
	
	_pKernel = new MeshKernelLinkedTriangles;
	clear();
	add_mesh(m);
	set_color(m.get_color());
	return *this;
}

void Mesh::add_mesh(const Mesh& f)
{
	// todo manage face ids

	//todo manage color by vertex
	_iColor = f.get_color();

	int iNbVertices = nb_vertices();
	for (int i = 0; i < f.nb_vertices(); i++)
	{
		Point3 v;
		f.get_vertex(i, v);
		_pKernel->add_vertex(v);
	}

	int iVertex1, iVertex2, iVertex3;
	for (int i = 0; i < f.nb_triangles(); i++)
	{
		if (f.is_triangle_unlinked(i))
			continue;

		f.get_triangle_vertices(i, iVertex1, iVertex2, iVertex3);
		add_triangle(iVertex1 + iNbVertices, iVertex2 + iNbVertices, iVertex3 + iNbVertices);
	}
}

int Mesh::add_triangle(int iVertex1, int iVertex2, int iVertex3)
{
	assert(iVertex1 >= 0);
	assert(iVertex2 >= 0);
	assert(iVertex3 >= 0);

	assert(iVertex1 < _pKernel->nb_vertices());
	assert(iVertex2 < _pKernel->nb_vertices());
	assert(iVertex3 < _pKernel->nb_vertices());

	return _pKernel->add_triangle(iVertex1,iVertex2,iVertex3);
}

int Mesh::add_triangle(const Point3& p1,const  Point3& p2,const Point3& p3)
{
	int t1 = _pKernel->add_vertex(p1);
	int t2 = _pKernel->add_vertex(p2);
	int t3 = _pKernel->add_vertex(p3);

	return add_triangle(t1, t2, t3);
}

void Mesh::add_quad(int iVertex1, int iVertex2, int iVertex3, int iVertex4)
{
	assert(iVertex1 >= 0);
	assert(iVertex2 >= 0);
	assert(iVertex3 >= 0);
	assert(iVertex4 >= 0);

	assert(iVertex1 < _pKernel->nb_vertices());
	assert(iVertex2 < _pKernel->nb_vertices());
	assert(iVertex3 < _pKernel->nb_vertices());
	assert(iVertex4 < _pKernel->nb_vertices());

	_pKernel->add_triangle(iVertex1, iVertex2, iVertex3);
	_pKernel->add_triangle(iVertex3, iVertex4, iVertex1);
}

void Mesh::add_quad(const Point3& p1, const  Point3& p2, const Point3& p3, const Point3& p4,bool bOptimSurface)
{
	int i1 = _pKernel->add_vertex(p1);
	int i2 = _pKernel->add_vertex(p2);
	int i3 = _pKernel->add_vertex(p3);
	int i4 = _pKernel->add_vertex(p4);

	if (!bOptimSurface)
		add_quad(i1, i2, i3, i4);
	else
	{
		//find the quad that as the smallest surface
		Triangle3 t1(p1, p2, p3);
		Triangle3 t2(p3, p4, p1);
		double surfaceQuad1 = t1.surface() + t2.surface();

		Triangle3 t3(p2, p3, p4);
		Triangle3 t4(p4, p1, p2);
		double surfaceQuad2 = t3.surface() + t4.surface();

		if (surfaceQuad1 <= surfaceQuad2)
			add_quad(i1, i2, i3, i4);
		else
			add_quad(i2, i3, i4, i1);
	}
}

void Mesh::add_pentagon(int iVertex1, int iVertex2, int iVertex3, int iVertex4, int iVertex5)
{
	assert(iVertex1 >= 0);
	assert(iVertex2 >= 0);
	assert(iVertex3 >= 0);
	assert(iVertex4 >= 0);
	assert(iVertex5 >= 0);

	assert(iVertex1 < _pKernel->nb_vertices());
	assert(iVertex2 < _pKernel->nb_vertices());
	assert(iVertex3 < _pKernel->nb_vertices());
	assert(iVertex4 < _pKernel->nb_vertices());
	assert(iVertex5 < _pKernel->nb_vertices());

	_pKernel->add_triangle(iVertex1, iVertex2, iVertex3);
	_pKernel->add_triangle(iVertex3, iVertex4, iVertex5);
	_pKernel->add_triangle(iVertex1, iVertex3, iVertex5);
}

void Mesh::split_triangle_with_vertex(int iTriangle, const Point3& p)
{
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	int iVertexP1 = add_vertex(p);
	split_triangle_with_vertex(iTriangle, iVertexP1);
}

void Mesh::split_triangle_with_vertex(int iTriangle, int iVertex)
{ 
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	assert(iVertex >= 0);
	assert(iVertex < _pKernel->nb_vertices());

	//replace triangle with 3 triangles built with iVertex
	int iV1, iV2, iV3;
	_pKernel->get_triangle_vertices(iTriangle, iV1, iV2, iV3);
	_pKernel->unlink_triangle(iTriangle);

	_pKernel->add_triangle(iV1, iV2, iVertex);
	_pKernel->add_triangle(iV2, iV3, iVertex);
	_pKernel->add_triangle(iV3, iV1, iVertex);
}
//////////////////////////////////////////////////////////////////////////////////
// split edge at new vertex
void Mesh::split_edge_with_vertex(int iTriangle1, int iTriangle2, int iVertex1, int iVertex2, int iVertexSplit) // 4 new triangles added at the end
{
	assert(iTriangle1 != -1);
	assert(iTriangle2 != -1);

	//split triangle1
	int iV1, iV2, iV3;
	_pKernel->get_triangle_vertices(iTriangle1, iV1, iV2, iV3);
	if (iV1 != iVertex1)
		rotate(iV1, iV2, iV3);
	if (iV1 != iVertex1)
		rotate(iV1, iV2, iV3);
	assert(iV1 == iVertex1);

	_pKernel->add_triangle(iV1, iVertexSplit, iV3);
	_pKernel->add_triangle(iVertexSplit, iV2, iV3);
	_pKernel->unlink_triangle(iTriangle1);

	//split triangle2
	_pKernel->get_triangle_vertices(iTriangle2, iV1, iV2, iV3);
	if (iV2 != iVertex2)
		rotate(iV1, iV2, iV3);
	if (iV2 != iVertex2)
		rotate(iV1, iV2, iV3);
	assert(iV2 == iVertex2);

	_pKernel->add_triangle(iV2, iVertexSplit, iV3);
	_pKernel->add_triangle(iVertexSplit, iV1, iV3);
	_pKernel->unlink_triangle(iTriangle2);
}
//////////////////////////////////////////////////////////////////////////////////
void Mesh::flip_triangle(int iTriangle)
{ 	
	assert(iTriangle >= 0);
	assert(iTriangle < _pKernel->nb_triangles());

	int iV1, iV2, iV3;
	_pKernel->get_triangle_vertices(iTriangle, iV1, iV2, iV3);
	_pKernel->unlink_triangle(iTriangle);
	_pKernel->add_triangle(iV1, iV3, iV2);
}
//////////////////////////////////////////////////////////////////////////////////
void Mesh::clear()
{
	_pKernel->clear();
}

bool Mesh::empty() const
{
	return _pKernel->nb_triangles()==0;
}

int Mesh::add_vertex(const Point3& vertex)
{
	return _pKernel->add_vertex(vertex);
}

int Mesh::add_vertex(double a, double b, double c)
{
	return _pKernel->add_vertex(Point3(a, b, c));
}

void Mesh::set_vertex(int iVertex, const Point3& vertex)
{
	assert(iVertex >= 0);
	assert(iVertex < _pKernel->nb_vertices());

	_pKernel->set_vertex(iVertex,vertex);
}

void Mesh::get_vertex(int iVertex, Point3& vertex) const
{
	assert(iVertex >= 0);
	assert(iVertex < _pKernel->nb_vertices());

	_pKernel->get_vertex(iVertex,vertex);
}

int Mesh::nb_vertices() const
{
	return _pKernel->nb_vertices();
}

void Mesh::apply_transform(const Transform& t)
{
	for (int i = 0; i < nb_vertices(); i++)
	{
		// todo optimize
		Point3 pv;
		get_vertex(i, pv);
		t.apply(pv);
		set_vertex(i, pv);
	}
}
//////////////////////////////////////////////////////////////////////////////////
void Mesh::split_triangle(int iTriangle, const Triangle3 & tSplitter)
{
	Triangle3 tA;
	get_triangle(iTriangle, tA);

	if (!tA.intersect_with(tSplitter))
		return;

	const double EPS = 1.e-8;
	vector<Point3> vIntersections;

	auto add_unique_intersection = [&](const Point3& p)
	{
		for (size_t i = 0; i < vIntersections.size(); ++i)
		{
			if ((vIntersections[i] - p).norm() <= EPS)
				return;
		}
		vIntersections.push_back(p);
	};

	auto collect_segment_triangle = [&](const Triangle3& tRef, const Segment3& edge)
	{
		Point3 p;
		if (tRef.intersect_with(edge, p))
			add_unique_intersection(p);
	};

	collect_segment_triangle(tA, Segment3(tSplitter.p1(), tSplitter.p2()));
	collect_segment_triangle(tA, Segment3(tSplitter.p2(), tSplitter.p3()));
	collect_segment_triangle(tA, Segment3(tSplitter.p3(), tSplitter.p1()));

	collect_segment_triangle(tSplitter, Segment3(tA.p1(), tA.p2()));
	collect_segment_triangle(tSplitter, Segment3(tA.p2(), tA.p3()));
	collect_segment_triangle(tSplitter, Segment3(tA.p3(), tA.p1()));

	if (vIntersections.empty())
		return;

	if (vIntersections.size() == 1)
	{
		split_triangle_with_vertex(iTriangle, vIntersections[0]);
		return;
	}

	Point3 P1 = vIntersections[0];
	Point3 P2 = vIntersections[1];
	double dMax = (P1 - P2).norm_square();
	for (size_t i = 0; i < vIntersections.size(); ++i)
	{
		for (size_t j = i + 1; j < vIntersections.size(); ++j)
		{
			double d = (vIntersections[i] - vIntersections[j]).norm_square();
			if (d > dMax)
			{
				dMax = d;
				P1 = vIntersections[i];
				P2 = vIntersections[j];
			}
		}
	}

	if (dMax <= EPS * EPS)
		return;

	int iVertice = add_vertex(P1);
	split_triangle_with_vertex(iTriangle, iVertice);

	int iVertice2 = add_vertex(P2);

	int iT1 = nb_triangles() - 3;
	int iT2 = nb_triangles() - 2;
	int iT3 = nb_triangles() - 1;

	Triangle3 tSplit;
	get_triangle(iT1, tSplit);
	if (tSplit.contains(P2))
	{
		split_triangle_with_vertex(iT1, iVertice2);
		return;
	}

	get_triangle(iT2, tSplit);
	if (tSplit.contains(P2))
	{
		split_triangle_with_vertex(iT2, iVertice2);
		return;
	}

	get_triangle(iT3, tSplit);
	if (tSplit.contains(P2))
		split_triangle_with_vertex(iT3, iVertice2);

}

