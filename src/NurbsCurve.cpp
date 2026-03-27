#include "NurbsCurve.h"

#include <cassert>
#include <algorithm>
#include <cmath>

#include "NurbsBasis.h"
#include "NurbsKnots.h"

///////////////////////////////////////////////////////////////////////////
NurbsCurve::NurbsCurve() :
	_degree(0), _iNbPoints(0)
{
}

NurbsCurve::~NurbsCurve()
{
}

void NurbsCurve::clear()
{
	_degree = 0;
	_knots.clear();
	_weights.clear();
	_points.clear();
	_iNbPoints = 0;
}

void NurbsCurve::set_degree(int degree)
{
	_degree = degree;
	_tempPoints.resize(degree + 1);
	_tempWeights.resize(degree + 1);
}

int NurbsCurve::degree() const
{
	return _degree;
}

void NurbsCurve::set_knots(const std::vector <double>& knots)
{
	_knots = knots;
	scale_knots(_knots);
}

void NurbsCurve::scale_knots(std::vector<double>& knots)
{
	NurbsKnots::normalize_to_01(knots);
}

void NurbsCurve::set_uniform()
{
	std::vector <double> newKnots;

	for (int i = 0; i <= _degree; i++)
		newKnots.push_back(0.);

	for (int i = 1; i < _iNbPoints - _degree; i++)
		newKnots.push_back(i);

	for (int i = 0; i <= _degree; i++)
		newKnots.push_back(_iNbPoints - _degree);

	set_knots(newKnots);
}

const std::vector<double>& NurbsCurve::knots() const
{
	return _knots;
}

void NurbsCurve::set_weights(const std::vector <double>& weights)
{
	_weights = weights;
}

const std::vector<double>& NurbsCurve::weights() const
{
	return _weights;
}

void NurbsCurve::set_equals_weights() //non rational
{
	_weights.assign(_iNbPoints, 1.);
}

void NurbsCurve::set_points(const std::vector <Point3>& points)
{
	_points = points;
	_iNbPoints = _points.size();
}

int NurbsCurve::nb_points() const
{
	return _points.size();
}

const std::vector<Point3>& NurbsCurve::points() const
{
	return _points;
}

std::vector<Point3>& NurbsCurve::points()
{
	return _points;
}

bool NurbsCurve::is_closed(double dTol) const
{
	if (_points.size() <= 1)
		return false;
	return (_points[0] - _points[_points.size() - 1]).norm_square() < (dTol * dTol);
}

int NurbsCurve::find_knot_span(const std::vector <double>& knots, double t)
{
	if (knots.size() < 2)
		return 0;

	if (t < 0.)
		t = 0.;

	if (t > 1.)
		t = 1.;

	//simple linear search for now
	for (int i = 0; i < knots.size() - 1; i++)
	{
		if (t < 1.)
		{
			if ((t >= knots[i]) && (t < knots[i + 1]))
				return i;
		}
		else
		{
			if ((t > knots[i]) && (t <= knots[i + 1]))
				return i;
		}
	}

	return (int)knots.size() - 2;

	/*
	//dichotomy
	int a = 0;
	int c = knots.size()-1;

	if (t < 0.)
		t = 0.;

	if (t > 1.)
		t = 1.;

	while (c-a > 1)
	{
		assert(a < c);

		int b = (a + c)/2;

		if (t >= knots[b])
			a = b;
		else
			c = b;

		assert(a >= 0);
		assert(a < knots.size());
		assert(c >= 0);
		assert(c < knots.size());

		assert(t >= knots[a]);
		assert(t <= knots[c]);
		assert(knots[a] < knots[c]);
	}

	assert(knots[a] < knots[a + 1]);
	assert(t < knots[a+1]);

	return a;
	*/
}

void NurbsCurve::insert_knot(double u)
{
	// no test for multiplicity for now
	int indexU = find_knot_span(_knots, u);

	std::vector<double> k = _knots;
	std::vector<Point3> p;
	std::vector<double> w;

	//reconstruct the curve
	for (int i = 0; i < _points.size() + 1; i++)
	{
		if (i <= indexU - _degree)
		{
			p.push_back(_points[i] * _weights[i]);
			w.push_back(_weights[i]);
			continue;
		}

		if (i >= indexU + 1)
		{
			p.push_back(_points[i - 1] * _weights[i - 1]);
			w.push_back(_weights[i - 1]);
			continue;
		}

		double dAlpha = (u - _knots[i]) / (_knots[i + _degree] - _knots[i]);
		p.push_back(_points[i] * _weights[i] * dAlpha + _points[i - 1] * _weights[i - 1] * (1. - dAlpha));
		w.push_back(_weights[i] * dAlpha + _weights[i - 1] * (1. - dAlpha));
	}

	//Homogenous to 3D coords
	for (int i = 0; i < p.size(); i++)
		p[i] /= w[i];

	k.insert(k.begin() + indexU + 1, u);
	set_knots(k);
	set_weights(w);
	set_points(p);
}

void NurbsCurve::evaluate(double u, Point3& p) const
{
	if (_iNbPoints == 0)
	{
		p = Point3();
		return;
	}

	//todo optimize all:
	assert(_iNbPoints == _weights.size());
	assert(_iNbPoints == _knots.size() - _degree - 1); //todo, always true ?

	int knotIndex = find_knot_span(_knots, u);

	for (int j = 0; j < _degree + 1; j++)
	{
		int idx = j + knotIndex - _degree;
		assert(idx >= 0);
		assert(idx < _iNbPoints);

		double w = _weights[idx];
		_tempWeights[j] = w;
		_tempPoints[j] = _points[idx] * w;
	}

	for (int r = 1; r < _degree + 1; r++)
		for (int j = _degree; j > r - 1; j--)
		{
			double denom = _knots[j + 1 + knotIndex - r] - _knots[j + knotIndex - _degree];
			double alpha = 0.;
			if (denom != 0.)
				alpha = (u - _knots[j + knotIndex - _degree]) / denom;
			_tempPoints[j] = _tempPoints[j - 1] * (1. - alpha) + _tempPoints[j] * alpha;
			_tempWeights[j] = _tempWeights[j - 1] * (1. - alpha) + _tempWeights[j] * alpha;
		}

	if (_tempWeights[_degree] == 0.)
		p = _tempPoints[_degree];
	else
		p = _tempPoints[_degree] / _tempWeights[_degree];
}

void NurbsCurve::evaluate_derivatives(double u, Point3& d1, Point3& d2) const
{
	d1 = Point3();
	d2 = Point3();

	if (_iNbPoints == 0 || _degree < 1)
		return;

	assert(_iNbPoints == (int)_weights.size());
	assert(_iNbPoints == (int)_knots.size() - _degree - 1);

	const int span = find_knot_span(_knots, u);
	double ders[3][32];
	NurbsBasis::basis_function_derivatives(_degree, _knots, span, u, ders);

	Point3 A0, A1, A2;
	double W0 = 0.;
	double W1 = 0.;
	double W2 = 0.;

	for (int j = 0; j <= _degree; ++j)
	{
		const int idx = span - _degree + j;
		if (idx < 0 || idx >= _iNbPoints)
			continue;

		const double w = _weights[idx];
		const Point3 Pw = _points[idx] * w;

		A0 += Pw * ders[0][j];
		A1 += Pw * ders[1][j];
		A2 += Pw * ders[2][j];

		W0 += w * ders[0][j];
		W1 += w * ders[1][j];
		W2 += w * ders[2][j];
	}

	const double wAbs = std::fabs(W0);
	if (wAbs < 1.e-14)
		return;

	const Point3 C = A0 / W0;
	const double W0sq = W0 * W0;
	const double W0cb = W0sq * W0;

	d1 = (A1 * W0 - A0 * W1) / W0sq;
	d2 = (A2 * W0sq - A0 * (W0 * W2) - A1 * (2. * W0 * W1) + A0 * (2. * W1 * W1)) / W0cb;
	(void)C;
}

bool NurbsCurve::tangent(double u, Point3& t) const
{
	Point3 d1, d2;
	evaluate_derivatives(u, d1, d2);
	const double n2 = d1.norm_square();
	if (n2 < 1.e-20)
	{
		t = Point3();
		return false;
	}
	t = d1 / std::sqrt(n2);
	return true;
}

bool NurbsCurve::normal(double u, Point3& n) const
{
	Point3 d1, d2;
	evaluate_derivatives(u, d1, d2);

	const double d1n2 = d1.norm_square();
	if (d1n2 < 1.e-20)
	{
		n = Point3();
		return false;
	}

	Point3 t = d1 / std::sqrt(d1n2);
	Point3 a = d2 - t * d2.dot_product(t);
	const double an2 = a.norm_square();
	if (an2 < 1.e-20)
	{
		n = Point3();
		return false;
	}

	n = a / std::sqrt(an2);
	return true;
}

double NurbsCurve::curvature(double u) const
{
	Point3 d1, d2;
	evaluate_derivatives(u, d1, d2);

	const double d1n = d1.norm();
	if (d1n < 1.e-20)
		return 0.;

	const Point3 c = d1.cross_product(d2);
	return c.norm() / (d1n * d1n * d1n);
}

///////////////////////////////////////////////////////////////////////////
void NurbsCurve::to_polyline(std::vector<Point3>& polyline) const
{
	polyline.clear();
	if (_points.empty())
		return;

	int iSamples = (int)_points.size() * 10;
	if (iSamples < 1)
		iSamples = 1;

	Point3 p;
	for (int i = 0; i <= iSamples; i++)
	{
		double t = (double)i / iSamples;
		evaluate(t, p);
		polyline.push_back(p);
	}
}
///////////////////////////////////////////////////////////////////////////
bool NurbsCurve::degree_elevation()
{
	if (_degree >= 3)
		return false; //only deg <=3 are handled

	// from paper: DIRECT DEGREE ELEVATION OF NURBS CURVES Kestutis Jankauskas, Dalius Rubliauskas

	//Procedure ElevateKnots
	std::vector<double> elevatedKnots, knotMultiplicity;
	int nbKnots = _knots.size();
	int lMult = 1;
	for (int k = 0; k < nbKnots - 1; k++)
	{
		elevatedKnots.push_back(_knots[k]);

		if (_knots[k] != _knots[k + 1])
		{
			elevatedKnots.push_back(_knots[k]);
			knotMultiplicity.push_back(lMult);
			lMult = 0;
		}
		lMult++;
	}
	knotMultiplicity.push_back(lMult);
	elevatedKnots.push_back(_knots[nbKnots - 1]);
	elevatedKnots.push_back(_knots[nbKnots - 1]);

	std::vector<Point3> newPoints;// (2 * _iNbPoints - 1);
	std::vector<double> newWeights;// (2 * _iNbPoints - 1);

	if (_degree == 1)
	{
		for (int i = 0; i <= _iNbPoints - 2; i++)
		{
			newPoints.push_back(_points[i] * _weights[i]);
			newWeights.push_back(_weights[i]);

			newPoints.push_back((_points[i] * _weights[i] + _points[i + 1] * _weights[i + 1]) / 2.);
			newWeights.push_back((_weights[i] + _weights[i + 1]) / 2.);
		}
		newPoints.push_back(_points[_iNbPoints - 1] * _weights[_iNbPoints - 1]);
		newWeights.push_back(_weights[_iNbPoints - 1]);
	}
	else if (_degree == 2)
	{
		//Procedure DDEQuadratic
		// in: U[m] � knot vector ==_knots
		// in: P[n] � control points ==_points
		// in: S[ns] � knot multiplicity vector ==knotMultiplicity
		// out: Q[n+ns-1] � elevated control points ==newPoints

		int k = 2;
		double b1 = 1. / 3.;
		double b2 = 2. / 3.;

		newPoints.push_back(_points[0] * _weights[0]);
		newWeights.push_back(_weights[0]);

		for (int l = 1; l < knotMultiplicity.size(); l++) // was //for (int l = 1 ; l<= nl � 1 ; l++)
		{
			if (knotMultiplicity[l - 1] > 1)
			{
				newPoints.push_back(_points[k - 1] * _weights[k - 1] * b2 + _points[k - 2] * _weights[k - 2] * b1);
				newWeights.push_back(_weights[k - 1] * b2 + _weights[k - 2] * b1);
			}
			if (knotMultiplicity[l] > 1)
			{
				newPoints.push_back(_points[k] * _weights[k] * b1 + _points[k - 1] * _weights[k - 1] * b2);
				newWeights.push_back(_weights[k] * b1 + _weights[k - 1] * b2);

				newPoints.push_back(_points[k] * _weights[k]);
				newWeights.push_back(_weights[k]);
			}
			else
			{
				double B = (_knots[k + 1] - _knots[k]) / (3 * (_knots[k + 2] - _knots[k]));

				newPoints.push_back(_points[k] * _weights[k] * B + _points[k - 1] * _weights[k - 1] * (1 - B));
				newWeights.push_back(_weights[k] * B + _weights[k - 1] * (1 - B));

				newPoints.push_back(_points[k] * _weights[k] * (B + b2) + _points[k - 1] * _weights[k - 1] * (b1 - B));
				newWeights.push_back(_weights[k] * (B + b2) + _weights[k - 1] * (b1 - B));
			}
			k += knotMultiplicity[l];
		}

		//	newPoints[2 * _iNbPoints - 2] = _points[_iNbPoints - 1] * _weights[_iNbPoints - 1];
		//	newWeights[2 * _iNbPoints - 2] = _weights[_iNbPoints - 1];
	}

	//homogenous to 3d points
	for (int i = 0; i < newPoints.size(); i++)
		newPoints[i] /= newWeights[i];

	set_knots(elevatedKnots);
	set_points(newPoints);
	set_weights(newWeights);
	set_degree(_degree + 1);
	return true;
}
///////////////////////////////////////////////////////////////////////////
namespace NurbsCurveUtil
{
	void to_polyline_adaptative_res(const NurbsCurve& curve, std::vector<Point3>& polyline, double dTol)
	{
		polyline.clear();
		if (curve.nb_points() == 0)
			return;

		std::vector<Point3> stack;
		stack.push_back(curve.points()[0]);
		stack.push_back(curve.points()[curve.nb_points() - 1]);

		while (!stack.empty())
		{
			Point3 p2 = stack.back();
			stack.pop_back();
			Point3 p1 = stack.back();
			stack.pop_back();

			double d = (p2 - p1).norm();
			if (d < dTol)
			{
				polyline.push_back(p1);
				continue;
			}

			double u = 0.5;
			Point3 pm;
			curve.evaluate(u, pm);

			stack.push_back(p2);
			stack.push_back(pm);
			stack.push_back(pm);
			stack.push_back(p1);
		}
	}


	bool has_valid_knot_vector(const std::vector<double>& knots, int degree, int nbPoints)
	{
		if (nbPoints <= 0)
			return false;

		if ((int)knots.size() != nbPoints + degree + 1)
			return false;

		for (int i = 0; i < (int)knots.size(); ++i)
		{
			if (!std::isfinite(knots[i]))
				return false;

			if (i > 0 && knots[i] < knots[i - 1])
				return false;
		}

		return true;
	}

	std::vector<double> build_clamped_uniform_knots(int nbPoints, int degree)
	{
		return NurbsKnots::build_clamped_uniform_knots(degree, nbPoints);
	}

	bool curves_are_compatible(const NurbsCurve& c1, const NurbsCurve& c2)
	{
		if (c1.nb_points() <= 0 || c2.nb_points() <= 0)
			return false;
		if (c1.nb_points() != c2.nb_points())
			return false;
		if ((int)c1.weights().size() != c1.nb_points() || (int)c2.weights().size() != c2.nb_points())
			return false;
		if (c1.degree() != c2.degree())
			return false;
		if (c1.knots().size() != c2.knots().size())
			return false;

		for (int i = 0; i < (int)c1.knots().size(); ++i)
			if (std::fabs(c1.knots()[i] - c2.knots()[i]) > 1.e-12)
				return false;

		return true;
	}

	void reverse_curve(NurbsCurve& c)
	{
		std::vector<Point3> points = c.points();
		std::vector<double> weights = c.weights();
		std::vector<double> knots = c.knots();

		std::reverse(points.begin(), points.end());
		std::reverse(weights.begin(), weights.end());

		if (!knots.empty())
		{
			const double kMin = knots.front();
			const double kMax = knots.back();
			std::vector<double> newKnots(knots.size(), 0.);
			for (int i = 0; i < (int)knots.size(); ++i)
				newKnots[i] = kMin + kMax - knots[(int)knots.size() - 1 - i];
			knots.swap(newKnots);
		}

		c.set_points(points);
		c.set_weights(weights);
		c.set_knots(knots);
	}

	bool align_curve_orientation(NurbsCurve& c1, NurbsCurve& c2)
	{
		if (c1.nb_points() != c2.nb_points() || c1.nb_points() <= 0)
			return false;

		const std::vector<Point3>& p1 = c1.points();
		const std::vector<Point3>& p2 = c2.points();
		const double dForward = p1.front().distance_square(p2.front()) + p1.back().distance_square(p2.back());
		const double dReverse = p1.front().distance_square(p2.back()) + p1.back().distance_square(p2.front());

		if (dReverse < dForward)
			NurbsCurveUtil::reverse_curve(c2);

		return true;
	}


	std::vector<double> build_open_uniform_knots(int degree, int nCtrl)
	{
		return NurbsKnots::build_open_uniform_knots(degree, nCtrl);
	}
}



