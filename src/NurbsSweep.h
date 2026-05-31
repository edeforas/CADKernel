#ifndef NurbsSweep_
#define NurbsSweep_

class NurbsCurve;
class NurbsSurface;
class NurbsSolid;

namespace NurbsSweep
{
	bool sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface, bool perpendicular = false);
	
	bool sweep_solid(const NurbsCurve& profile, const NurbsCurve& path, NurbsSolid& solid, bool perpendicular = true);
}

#endif
