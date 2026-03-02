#ifndef NurbsFitting_
#define NurbsFitting_

#include <vector>

class Point3;
class NurbsCurve;
class NurbsSurface;

namespace NurbsFitting
{
	bool fit_curve_least_squares(
		const std::vector<Point3>& samples,
		int degree,
		int iNbCtrl,
		NurbsCurve& outCurve,
		double dRegularization = 1.e-8);

	bool fit_surface_least_squares(
		const std::vector<Point3>& samples,
		int iNbSamplesU,
		int iNbSamplesV,
		int degreeU,
		int degreeV,
		int iNbCtrlU,
		int iNbCtrlV,
		NurbsSurface& outSurface,
		double dRegularization = 1.e-8);
}

#endif
