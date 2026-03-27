#include "NurbsOffset.h"

#include "NurbsSurface.h"

namespace NurbsOffset
{
	void offset_surface(const NurbsSurface& input, double distance, NurbsSurface& output)
	{
		output = input; // copy all properties

		const std::vector<Point3>& points = input.points();
		std::vector<Point3>& out_points = output.points();
		int nbU = input.nb_points_u();
		int nbV = input.nb_points_v();

		const std::vector<double>& knotsU = input.knots_u();
		const std::vector<double>& knotsV = input.knots_v();

		// For each control point, compute u,v and normal
		for (int iV = 0; iV < nbV; ++iV)
		{
			for (int iU = 0; iU < nbU; ++iU)
			{
				int idx = iV * nbU + iU;
				// Approximate u,v for control point
				// Assuming uniform distribution
				double u = (double)iU / (nbU - 1);
				double v = (double)iV / (nbV - 1);
				// Clamp to [0,1] and scale to knot range

				if (!knotsU.empty() && !knotsV.empty())
				{
					double uMin = knotsU.front();
					double uMax = knotsU.back();
					double vMin = knotsV.front();
					double vMax = knotsV.back();
					u = uMin + u * (uMax - uMin);
					v = vMin + v * (vMax - vMin);
				}

				Point3 normal;
				if (input.normal(u, v, normal))
				{
					out_points[idx] = points[idx] + normal * distance;
				}
				else
				{
					// If normal fails, offset in Z direction
					out_points[idx] = points[idx] + Point3(0, 0, distance);
				}
			}
		}
	}
}