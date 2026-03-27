#ifndef _NurbsCurveJoin_
#define _NurbsCurveJoin_

#include "NurbsCurve.h"

class NurbsCurveJoin
{
public:
	enum Continuity
	{
		G0 = 0,
		G1 = 1,
		G2 = 2
	};

	static bool create_connector(const NurbsCurve& c1, const NurbsCurve& c2, Continuity continuity, NurbsCurve& out);
	static bool create_chamfer(const NurbsCurve& c1, const NurbsCurve& c2, double chamferLength, NurbsCurve& out);
	static bool create_fillet(const NurbsCurve& c1, const NurbsCurve& c2, double radius, NurbsCurve& out);
};

#endif
