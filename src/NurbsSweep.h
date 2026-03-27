#ifndef NurbsSweep_
#define NurbsSweep_

class NurbsCurve;
class NurbsSurface;

namespace NurbsSweep
{
	bool sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface, bool perpendicular = false);
}

#endif
