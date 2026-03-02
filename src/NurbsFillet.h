#ifndef NurbsFillet_
#define NurbsFillet_

class NurbsCurve;
class NurbsSurface;
class NurbsSolid;

class NurbsFillet
{
public:
	enum SurfaceEdge
	{
		EdgeUMin,
		EdgeUMax,
		EdgeVMin,
		EdgeVMax
	};

	NurbsFillet();
	virtual ~NurbsFillet();

	bool create_chamfer(const NurbsCurve& c1, const NurbsCurve& c2, double dOffset1, double dOffset2, NurbsSurface& out) const;
	bool create_fillet(const NurbsCurve& c1, const NurbsCurve& c2, double dRadius, NurbsSurface& out) const;

	bool create_chamfer_between_surfaces(const NurbsSurface& s1, SurfaceEdge e1, const NurbsSurface& s2, SurfaceEdge e2, double dOffset1, double dOffset2, NurbsSurface& out) const;
	bool create_fillet_between_surfaces(const NurbsSurface& s1, SurfaceEdge e1, const NurbsSurface& s2, SurfaceEdge e2, double dRadius, NurbsSurface& out) const;
	bool find_shared_edge_pair(const NurbsSurface& s1, const NurbsSurface& s2, SurfaceEdge& e1, SurfaceEdge& e2, double dTol = 1.e-4) const;

	bool create_chamfer_on_solid(const NurbsSolid& solid, int iSurfaceA, SurfaceEdge eA, int iSurfaceB, SurfaceEdge eB, double dOffset1, double dOffset2, NurbsSurface& out) const;
	bool create_fillet_on_solid(const NurbsSolid& solid, int iSurfaceA, SurfaceEdge eA, int iSurfaceB, SurfaceEdge eB, double dRadius, NurbsSurface& out) const;
	bool create_chamfer_on_solid_auto(const NurbsSolid& solid, int iSurfaceA, int iSurfaceB, double dOffset1, double dOffset2, NurbsSurface& out) const;
	bool create_fillet_on_solid_auto(const NurbsSolid& solid, int iSurfaceA, int iSurfaceB, double dRadius, NurbsSurface& out) const;
};

#endif
