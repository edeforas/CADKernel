#ifndef NurbsSurface_
#define NurbsSurface_

#include <vector>

#include "Geometry.h"
#include "NurbsCurve.h"
class Transform;
class NurbsTrimmedSurface;

///////////////////////////////////////////////////////////////////////////
class NurbsSurface
{
public:
	NurbsSurface();
	virtual ~NurbsSurface();
	NurbsSurface& operator=(const NurbsSurface& other);

	void clear();

	void set_degree(int degreeU, int degreeV);
	int degree_u() const;
	int degree_v() const;

	void set_knots_u(const std::vector <double>& knots);
	void set_knots_v(const std::vector <double>& knots);
	void set_uniform_u();
	void set_uniform_v();
	const std::vector<double>& knots_u() const;
	const std::vector<double>& knots_v() const;

	void set_weights(const std::vector <double>& weights);
	void set_equals_weights(); //non rational
	const std::vector<double>& weights() const;
	std::vector<double>& weights();

	void set_closed_u(bool bClosedU);
	void set_closed_v(bool bClosedV);

	void set_points(const std::vector <Point3>& points, int iNbPointsU, int iNbPointsV);
	const std::vector<Point3>& points() const;
	std::vector<Point3>& points();
	int nb_points_u() const;
	int nb_points_v() const;
	int nb_points() const;

	void apply_transform(const Transform& t);

	bool is_closed_u() const;
	bool is_closed_v() const;

	void insert_knot_u(double u);
	void insert_knot_v(double v);
	void insert_knot_uv(double u, double v);

	NurbsCurve edge_u0() const;
	NurbsCurve edge_u1() const;
	NurbsCurve edge_v0() const;
	NurbsCurve edge_v1() const;

	bool degree_elevation_u();
	bool degree_elevation_v();
	void reverse_u();
	NurbsSurface reversed_u() const;

	void extend_u(double distance, bool extend_start = true);
	void extend_v(double distance, bool extend_start = true);
	NurbsSurface extended_u(double distance, bool extend_start = true) const;
	NurbsSurface extended_v(double distance, bool extend_start = true) const;

	void evaluate(double u, double v, Point3& p) const;
	void evaluate_clamped(double u, double v, Point3& p) const;
	void evaluate_derivatives(double u, double v, Point3& du, Point3& dv, Point3& duu, Point3& duv, Point3& dvv) const;
	void evaluate_partials(double u, double v, Point3& du, Point3& dv) const;
	bool normal(double u, double v, Point3& n) const;
	bool curvature(double u, double v, double& gaussian, double& mean) const;
	void project_point_on_surface(const Point3& target, double& u, double& v, Point3& projected) const;

	virtual bool is_trimmed() const;
	virtual NurbsTrimmedSurface* trimming();
	virtual const NurbsTrimmedSurface* trimming() const;

private:
	static int find_knot_span(const std::vector <double>& knots, double u);
	static void scale_knots(std::vector<double>& knots);

	void create_u_curves(std::vector<NurbsCurve>& vu) const;
	void create_v_curves(std::vector<NurbsCurve>& vv) const;

	void from_u_curves(const std::vector<NurbsCurve>& vu); //reuse V knots and degreee
	void from_v_curves(const std::vector<NurbsCurve>& vv); //reuse U knots and degreee

	int _degreeU, _degreeV;
	int _iNbPointsU, _iNbPointsV;
	bool _bClosedU, _bClosedV;

	std::vector <double> _knotsU, _knotsV;
	std::vector <double> _weights;
	std::vector <Point3> _points;

	mutable std::vector<Point3> _tempPointsU;
	mutable std::vector<Point3> _tempPointsV;
	mutable std::vector<double> _tempWeightsU;
	mutable std::vector<double> _tempWeightsV;
};
///////////////////////////////////////////////////////////////////////////

#endif