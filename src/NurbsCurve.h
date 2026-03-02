#ifndef _NurbsCurve_
#define _NurbsCurve_

#include <vector>

#include "Geometry.h"

///////////////////////////////////////////////////////////////////////////
class NurbsCurve
{
public:
	NurbsCurve();
	virtual ~NurbsCurve();

	void clear();

	void set_degree(int degree);
	int degree() const;

	void set_knots(const std::vector <double>& knots);
	const std::vector<double>& knots() const;
	void set_uniform();

	void set_weights(const std::vector <double>& weights);
	const std::vector<double>& weights() const;
	void set_equals_weights(); //non rational

	void set_points(const std::vector <Point3>& points);
	int nb_points() const;

	bool is_closed(double dTol = 1.e-6) const;

	void insert_knot(double u);
	bool degree_elevation();

	const std::vector<Point3>& points() const;
	std::vector<Point3>& points();
	void evaluate(double u, Point3& p) const;
	void evaluate_derivatives(double u, Point3& d1, Point3& d2) const;
	bool tangent(double u, Point3& t) const;
	bool normal(double u, Point3& n) const;
	double curvature(double u) const;

	void to_polyline(std::vector<Point3>& polyline) const;

	static int find_knot_span(const std::vector <double>& knots, double t);

private:
	static void scale_knots(std::vector<double>& knots);

	int _degree;
	std::vector <double> _knots;
	std::vector <double> _weights;
	std::vector <Point3> _points;
	int _iNbPoints;

	mutable std::vector<Point3> _tempPoints;
	mutable std::vector<double> _tempWeights;
};

namespace NurbsCurveUtil
{
	void to_polyline_adaptative_res(const NurbsCurve& curve, std::vector<Point3>& polyline, double dTol = 1.e-3);
	bool has_valid_knot_vector(const std::vector<double>& knots, int degree, int nbPoints);
	std::vector<double> build_clamped_uniform_knots(int nbPoints, int degree);
	bool curves_are_compatible(const NurbsCurve& c1, const NurbsCurve& c2);
	void reverse_curve(NurbsCurve& c);
	bool align_curve_orientation(NurbsCurve& c1, NurbsCurve& c2);
	std::vector<double> build_open_uniform_knots(int degree, int nCtrl);
}

///////////////////////////////////////////////////////////////////////////

#endif