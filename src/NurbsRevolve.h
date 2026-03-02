#ifndef NurbsRevolve_
#define NurbsRevolve_

#include "Geometry.h"

class NurbsCurve; 
class NurbsSurface;

///////////////////////////////////////////////////////////////////////////
class NurbsRevolve
{
public:
	NurbsRevolve();
	virtual ~NurbsRevolve();
	bool revolve(const NurbsCurve& nc, NurbsSurface& ns) const;
	bool revolve(const NurbsCurve& nc, double dAngleRad, NurbsSurface& ns) const;
	bool revolve_deg(const NurbsCurve& nc, double dAngleDeg, NurbsSurface& ns) const;
};
///////////////////////////////////////////////////////////////////////////

#endif