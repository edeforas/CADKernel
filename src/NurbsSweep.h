#ifndef NurbsSweep_
#define NurbsSweep_

class NurbsCurve;
class NurbsSurface;
class NurbsSolid;

namespace NurbsSweep
{
	enum class EndpointClosure
	{
		None = 0,      // No closure
		Disk = 1,      // Close with flat disk using the profile curve
		HalfSphere = 2 // Close with half-sphere
	};

	bool sweep(const NurbsCurve& profile, const NurbsCurve& path, NurbsSurface& surface, bool perpendicular = false, EndpointClosure startCap = EndpointClosure::None, EndpointClosure endCap = EndpointClosure::None);
	
	bool sweep_solid(const NurbsCurve& profile, const NurbsCurve& path, NurbsSolid& solid, bool perpendicular = false, EndpointClosure startCap = EndpointClosure::None, EndpointClosure endCap = EndpointClosure::None);
}

#endif
