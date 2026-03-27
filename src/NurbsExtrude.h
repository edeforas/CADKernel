#ifndef NurbsExtrude_
#define NurbsExtrude_

#include "Geometry.h"

class NurbsCurve; 
class NurbsSurface;
class NurbsSolid;

///////////////////////////////////////////////////////////////////////////
class NurbsExtrude
{
public:
	NurbsExtrude();
	virtual ~NurbsExtrude();
	bool extrude(const NurbsCurve& nc, const Point3& direction, NurbsSurface& ns) const;
	bool extrude_face(NurbsSolid& solid, int faceIndex, const Point3& direction) const;

};
///////////////////////////////////////////////////////////////////////////

#endif