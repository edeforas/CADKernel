#ifndef _NurbsBooleanFallback_
#define _NurbsBooleanFallback_

class NurbsSolid;
class Mesh;

///////////////////////////////////////////////////////////////////////////
class NurbsBooleanFallback
{
public:
	// Mesh fallback for overlap cases where NurbsBoolean cannot build a solid result yet.
	static void union_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments = 24);
	static void intersection_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments = 24);
	static void difference_mesh(const NurbsSolid& a, const NurbsSolid& b, Mesh& result, int iNbSegments = 24);
};
///////////////////////////////////////////////////////////////////////////

#endif