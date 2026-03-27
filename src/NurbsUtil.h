#ifndef NurbsUtil_
#define NurbsUtil_

#include <vector>

class Point3;
class NurbsCurve; 
class NurbsSurface;
class NurbsSolid;
class NurbsTrimmedSurface;
class BezierSurface;
class Mesh;

///////////////////////////////////////////////////////////////////////////
namespace NurbsUtil
{
	struct KnotAnalysis
	{
		std::vector<double> unique_knots;
		std::vector<int> multiplicities;
	};

	void create_curve_from_points(const std::vector<Point3>& points, int iDegree, NurbsCurve& n);
	void create_surface_from_z(const std::vector<double>& z, int iSizeX, int iSizeY, int iDegree, NurbsSurface& n);
	void create_solid_from_mesh(const Mesh & m, NurbsSolid& n);
	void to_mesh(const NurbsSurface& n, Mesh& m, int iNbSegments=5, bool bClearMesh=true);
	void to_mesh(const NurbsSolid& ns, Mesh& m, int iNbSegments=5);
	void to_mesh(const NurbsTrimmedSurface& ts, Mesh& m, int iNbSegments = 24, bool bClearMesh = true);
	void to_mesh(const std::vector<NurbsTrimmedSurface>& trimmedSurfaces, Mesh& m, int iNbSegments = 24);
	void to_controlpoints_mesh(const NurbsSurface& n, Mesh& m);
	void to_bezier_patches(const NurbsSurface& n, std::vector<BezierSurface>& patches);
	void solid_to_trimmed_surfaces(const NurbsSolid& src, std::vector<NurbsTrimmedSurface>& dst);

	double sanitize_weight(double value,double kEpsilonWeight = 1.e-12);
	std::vector<double> build_safe_weights(const std::vector<double>& weights, int expectedSize);
	std::vector<double> build_segmented_quadratic_knots(int nbSegments);
	std::vector<double> build_uniform_knots(int degree, int nbCtrlPoints);

	KnotAnalysis analyze_knots(const std::vector<double>& knots, int expectedDegree, int expectedCtrlPoints);



}
///////////////////////////////////////////////////////////////////////////

#endif