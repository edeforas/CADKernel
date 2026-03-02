#include "NurbsSurface.h"
#include "Transform.h"

#include <cassert>
#include <algorithm>
#include <cmath>

#include "NurbsCurve.h"

namespace
{
	void basis_function_derivatives(int degree, const std::vector<double>& knots, int span, double u, std::vector<std::vector<double>>& ders)
	{
		ders.assign(3, std::vector<double>(degree + 1, 0.));
		if (degree < 0)
			return;

		std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1, 0.));
		std::vector<double> left(degree + 1, 0.);
		std::vector<double> right(degree + 1, 0.);
		ndu[0][0] = 1.;

		for (int j = 1; j <= degree; ++j)
		{
			left[j] = u - knots[span + 1 - j];
			right[j] = knots[span + j] - u;
			double saved = 0.;
			for (int r = 0; r < j; ++r)
			{
				ndu[j][r] = right[r + 1] + left[j - r];
				double temp = 0.;
				if (std::fabs(ndu[j][r]) > 1.e-14)
					temp = ndu[r][j - 1] / ndu[j][r];
				ndu[r][j] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			ndu[j][j] = saved;
		}

		for (int j = 0; j <= degree; ++j)
			ders[0][j] = ndu[j][degree];

		std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1, 0.));
		for (int r = 0; r <= degree; ++r)
		{
			int s1 = 0;
			int s2 = 1;
			a[0][0] = 1.;

			for (int k = 1; k <= 2; ++k)
			{
				double d = 0.;
				int rk = r - k;
				int pk = degree - k;

				if (r >= k)
				{
					a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
					d = a[s2][0] * ndu[rk][pk];
				}

				int j1 = (rk >= -1) ? 1 : -rk;
				int j2 = (r - 1 <= pk) ? (k - 1) : (degree - r);
				for (int j = j1; j <= j2; ++j)
				{
					a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
					d += a[s2][j] * ndu[rk + j][pk];
				}

				if (r <= pk)
				{
					a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
					d += a[s2][k] * ndu[r][pk];
				}

				ders[k][r] = d;
				std::swap(s1, s2);
			}
		}

		int r = degree;
		for (int k = 1; k <= 2; ++k)
		{
			for (int j = 0; j <= degree; ++j)
				ders[k][j] *= r;
			r *= (degree - k);
		}
	}
}

///////////////////////////////////////////////////////////////////////////
NurbsSurface::NurbsSurface() :
	_degreeU(0), _degreeV(0), _iNbPointsU(0), _iNbPointsV(0), _bClosedU(false), _bClosedV(false)
{
}

NurbsSurface& NurbsSurface::operator=(const NurbsSurface& other)
{
	if (this == &other)
		return *this;

	_degreeU = other._degreeU;
	_degreeV = other._degreeV;
	_knotsU = other._knotsU;
	_knotsV = other._knotsV;
	_weights = other._weights;
	_points = other._points;
	_iNbPointsU = other._iNbPointsU;
	_iNbPointsV = other._iNbPointsV;
	_bClosedU = other._bClosedU;
	_bClosedV = other._bClosedV;

	return *this;
}

NurbsSurface::~NurbsSurface()
{
}

void NurbsSurface::clear()
{
	_degreeU = 0;
	_degreeV = 0;
	_knotsU.clear();
	_knotsV.clear();
	_weights.clear();
	_points.clear();
	_iNbPointsU = 0;
	_iNbPointsV = 0;
}

void NurbsSurface::set_degree(int degreeU, int degreeV)
{
	_degreeU = degreeU;
	_degreeV = degreeV;

	_tempPointsU.resize(degreeU + 1);
	_tempPointsV.resize(degreeV + 1);
	_tempWeightsU.resize(degreeU + 1);
	_tempWeightsV.resize(degreeV + 1);
}

int NurbsSurface::degree_u() const
{
	return _degreeU;
}

int NurbsSurface::degree_v() const
{
	return _degreeV;
}

void NurbsSurface::set_knots_u(const std::vector <double>& knots)
{
	_knotsU = knots;
	scale_knots(_knotsU);
}

void NurbsSurface::set_uniform_u()
{
	std::vector <double> knots;

	knots.clear();

	for (int i = 0; i <= _degreeU; i++)
		knots.push_back(0.);

	for (int i = 1; i < _iNbPointsU - _degreeU; i++)
		knots.push_back(i);

	for (int i = 0; i <= _degreeU; i++)
		knots.push_back(_iNbPointsU - _degreeU);

	set_knots_u(knots);
}

void NurbsSurface::set_knots_v(const std::vector <double>& knots)
{
	_knotsV = knots;
	scale_knots(_knotsV);
}

void NurbsSurface::set_uniform_v()
{
	std::vector <double> knots;

	for (int i = 0; i <= _degreeV; i++)
		knots.push_back(0.);

	for (int i = 1; i < _iNbPointsV - _degreeV; i++)
		knots.push_back(i);

	for (int i = 0; i <= _degreeV; i++)
		knots.push_back(_iNbPointsV - _degreeV);

	set_knots_v(knots);
}

const std::vector<double>& NurbsSurface::knots_u() const
{
	return _knotsU;
}
const std::vector<double>& NurbsSurface::knots_v() const
{
	return _knotsV;
}

void NurbsSurface::scale_knots(std::vector<double>& knots)
{
	if (knots.empty())
		return;

	//compute min and max
	double dMin = knots[0];
	double dMax = knots[0];
	for (int i = 1; i < knots.size(); i++)
	{
		double v = knots[i];
		dMin = std::min(dMin, v);
		dMax = std::max(dMax, v);
	}

	if (dMax == dMin)
	{
		for (int i = 0; i < knots.size(); i++)
			knots[i] = 0.;
		return;
	}

	//apply scale so that 0<= knots <= 1
	for (int i = 0; i < knots.size(); i++)
		knots[i] = (knots[i] - dMin) / (dMax - dMin);
}

void NurbsSurface::set_weights(const std::vector <double>& weights)
{
	_weights = weights;
}

void NurbsSurface::set_equals_weights() //non rational
{
	_weights.resize(_points.size(), 1.);
}

const std::vector<double>& NurbsSurface::weights() const
{
	return _weights;
}

std::vector<double>& NurbsSurface::weights()
{
	return _weights;
}

void NurbsSurface::set_points(const std::vector <Point3>& points, int iNbPointsU, int iNbPointsV)
{
	_points = points;
	_iNbPointsU = iNbPointsU;
	_iNbPointsV = iNbPointsV;
}

int NurbsSurface::nb_points_u() const
{
	return _iNbPointsU;
}

int NurbsSurface::nb_points_v() const
{
	return _iNbPointsV;
}

int NurbsSurface::nb_points() const
{
	return _points.size();
}

void NurbsSurface::apply_transform(const Transform& t)
{
	for (int i = 0; i < _points.size(); i++)
	{
		t.apply(_points[i]);
	}
}

void NurbsSurface::set_closed_u(bool bClosedU)
{
	_bClosedU = bClosedU;
}

void NurbsSurface::set_closed_v(bool bClosedV)
{
	_bClosedV = bClosedV;
}

bool NurbsSurface::is_closed_u() const
{
	return _bClosedU;
}

bool NurbsSurface::is_closed_v() const
{
	return _bClosedV;
}

const std::vector<Point3>& NurbsSurface::points() const
{
	return _points;
}

std::vector<Point3>& NurbsSurface::points()
{
	return _points;
}

int NurbsSurface::find_knot_span(const std::vector <double>& knots, double t)
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
}

void NurbsSurface::insert_knot_u(double u)
{
	std::vector<NurbsCurve> vu;
	create_u_curves(vu);

	for (int i = 0; i < vu.size(); i++)
		vu[i].insert_knot(u);

	from_u_curves(vu);
}

void NurbsSurface::insert_knot_v(double v)
{
	std::vector<NurbsCurve> vv;
	create_v_curves(vv);

	for (int i = 0; i < vv.size(); i++)
		vv[i].insert_knot(v);

	from_v_curves(vv);
}

void NurbsSurface::insert_knot_uv(double u, double v)
{
	insert_knot_u(u);
	insert_knot_v(v);
}

bool NurbsSurface::degree_elevation_u()
{
	if (_degreeU >= 3)
		return false; //only deg <=3 are handled

	std::vector<NurbsCurve> vu;
	create_u_curves(vu);

	for (int i = 0; i < vu.size(); i++)
		vu[i].degree_elevation();

	from_u_curves(vu);
	return true;
}

bool NurbsSurface::degree_elevation_v()
{
	if (_degreeV >= 3)
		return false; //only deg <=3 are handled

	std::vector<NurbsCurve> vv;
	create_v_curves(vv);

	for (int i = 0; i < vv.size(); i++)
		vv[i].degree_elevation();

	from_v_curves(vv);
	return true;
}

void NurbsSurface::create_u_curves(std::vector<NurbsCurve>& vu) const
{
	vu.clear();

	for (int v = 0; v < nb_points_v(); v++)
	{
		NurbsCurve n;
		n.set_degree(_degreeU);
		n.set_knots(_knotsU);

		std::vector<Point3> vp;
		std::vector<double> vw;
		for (int u = 0; u < nb_points_u(); u++)
		{
			vp.push_back(points()[u + nb_points_u() * v]);
			vw.push_back(weights()[u + nb_points_u() * v]);
		}

		n.set_points(vp);
		n.set_weights(vw);
		vu.push_back(n);
	}
}

void NurbsSurface::create_v_curves(std::vector<NurbsCurve>& vv) const
{
	vv.clear();

	for (int u = 0; u < nb_points_u(); u++)
	{
		NurbsCurve n;
		n.set_degree(_degreeV);
		n.set_knots(_knotsV);

		std::vector<Point3> vp;
		std::vector<double> vw;
		for (int v = 0; v < nb_points_v(); v++)
		{
			vp.push_back(points()[u + nb_points_u() * v]);
			vw.push_back(weights()[u + nb_points_u() * v]);
		}

		n.set_points(vp);
		n.set_weights(vw);
		vv.push_back(n);
	}
}

void NurbsSurface::from_u_curves(const std::vector<NurbsCurve>& vu) //reuse V knots and degree
{
	if (vu.empty())
		return;

	//update knots
	set_knots_u(vu[0].knots());

	//update_points
	std::vector<Point3> vp;
	std::vector<double> vw;
	for (int v = 0; v < vu.size(); v++)
	{
		vp.insert(vp.end(), vu[v].points().begin(), vu[v].points().end());
		vw.insert(vw.end(), vu[v].weights().begin(), vu[v].weights().end());
	}
	set_points(vp, vu[0].nb_points(), nb_points_v());
	set_weights(vw);
	set_degree(vu[0].degree(), degree_v());
}

void NurbsSurface::from_v_curves(const std::vector<NurbsCurve>& vv) //reuse U knots and degree
{
	if (vv.empty())
		return;

	//update knots
	set_knots_v(vv[0].knots());

	//update_points
	int nbPointsV = vv[0].nb_points();
	std::vector<Point3> vp(nb_points_u() * nbPointsV);
	std::vector<double> vw(nb_points_u() * nbPointsV);
	for (int u = 0; u < nb_points_u(); u++)
	{
		const NurbsCurve& nv = vv[u];

		for (int v = 0; v < nbPointsV; v++)
		{
			vp[u + v * nb_points_u()] = nv.points()[v];
			vw[u + v * nb_points_u()] = nv.weights()[v];
		}
	}
	set_points(vp, nb_points_u(), vv[0].nb_points());
	set_weights(vw);
	set_degree(degree_u(), vv[0].degree());
}

void NurbsSurface::evaluate(double u, double v, Point3& p) const
{
	if (_points.empty())
	{
		p = Point3();
		return;
	}

	//todo optimize all:
	assert(_points.size() == _weights.size());

	int knotIndexU = find_knot_span(_knotsU, u);
	int knotIndexV = find_knot_span(_knotsV, v);
	int iNbCtrlPointsU = nb_points_u();
	int iNbCtrlPointsV = nb_points_v();

	assert(iNbCtrlPointsU == _knotsU.size() - _degreeU - 1);
	assert(iNbCtrlPointsV == _knotsV.size() - _degreeV - 1);
	assert(iNbCtrlPointsU * iNbCtrlPointsV == _points.size());

	//tensor product
	for (int vi = 0; vi < _degreeV + 1; vi++)
	{
		int idxV = vi + knotIndexV - _degreeV;
		assert(idxV >= 0);
		assert(idxV < iNbCtrlPointsV);

		// evaluate on u direction
		for (int ui = 0; ui < _degreeU + 1; ui++)
		{
			int idxU = ui + knotIndexU - _degreeU;
			assert(idxU >= 0);
			assert(idxU < iNbCtrlPointsU);

			double w = _weights[idxU + iNbCtrlPointsU * idxV];
			_tempWeightsU[ui] = w;
			_tempPointsU[ui] = _points[idxU + iNbCtrlPointsU * idxV] * w;
		}

		for (int ru = 1; ru < _degreeU + 1; ru++)
			for (int ju = _degreeU; ju > ru - 1; ju--)
			{
				double denom = _knotsU[ju + 1 + knotIndexU - ru] - _knotsU[ju + knotIndexU - _degreeU];
				double alpha = 0.;
				if (denom != 0.)
					alpha = (u - _knotsU[ju + knotIndexU - _degreeU]) / denom;
				_tempPointsU[ju] = _tempPointsU[ju - 1] * (1. - alpha) + _tempPointsU[ju] * alpha;
				_tempWeightsU[ju] = _tempWeightsU[ju - 1] * (1. - alpha) + _tempWeightsU[ju] * alpha;
			}

		_tempPointsV[vi] = _tempPointsU[_degreeU];
		_tempWeightsV[vi] = _tempWeightsU[_degreeU];
	}

	//evaluate on v direction
	for (int rv = 1; rv < _degreeV + 1; rv++)
		for (int jv = _degreeV; jv > rv - 1; jv--)
		{
			double denom = _knotsV[jv + 1 + knotIndexV - rv] - _knotsV[jv + knotIndexV - _degreeV];
			double alpha = 0.;
			if (denom != 0.)
				alpha = (v - _knotsV[jv + knotIndexV - _degreeV]) / denom;
			_tempPointsV[jv] = _tempPointsV[jv - 1] * (1. - alpha) + _tempPointsV[jv] * alpha;
			_tempWeightsV[jv] = _tempWeightsV[jv - 1] * (1. - alpha) + _tempWeightsV[jv] * alpha;
		}

	if (_tempWeightsV[_degreeV] == 0.)
		p = _tempPointsV[_degreeV];
	else
		p = _tempPointsV[_degreeV] / _tempWeightsV[_degreeV];
}

void NurbsSurface::evaluate_derivatives(double u, double v, Point3& du, Point3& dv, Point3& duu, Point3& duv, Point3& dvv) const
{
	du = Point3();
	dv = Point3();
	duu = Point3();
	duv = Point3();
	dvv = Point3();

	if (_points.empty())
		return;

	assert(_iNbPointsU == (int)_knotsU.size() - _degreeU - 1);
	assert(_iNbPointsV == (int)_knotsV.size() - _degreeV - 1);
	assert(_iNbPointsU * _iNbPointsV == (int)_points.size());
	assert(_iNbPointsU * _iNbPointsV == (int)_weights.size());

	const int spanU = find_knot_span(_knotsU, u);
	const int spanV = find_knot_span(_knotsV, v);

	std::vector<std::vector<double>> Nu, Nv;
	basis_function_derivatives(_degreeU, _knotsU, spanU, u, Nu);
	basis_function_derivatives(_degreeV, _knotsV, spanV, v, Nv);

	Point3 A;
	Point3 Au;
	Point3 Av;
	Point3 Auu;
	Point3 Auv;
	Point3 Avv;
	double W = 0.;
	double Wu = 0.;
	double Wv = 0.;
	double Wuu = 0.;
	double Wuv = 0.;
	double Wvv = 0.;

	for (int j = 0; j <= _degreeV; ++j)
	{
		const int idxV = spanV - _degreeV + j;
		if (idxV < 0 || idxV >= _iNbPointsV)
			continue;

		for (int i = 0; i <= _degreeU; ++i)
		{
			const int idxU = spanU - _degreeU + i;
			if (idxU < 0 || idxU >= _iNbPointsU)
				continue;

			const int idx = idxU + _iNbPointsU * idxV;
			const double ww = _weights[idx];
			const Point3 Pw = _points[idx] * ww;

			const double b00 = Nu[0][i] * Nv[0][j];
			const double b10 = Nu[1][i] * Nv[0][j];
			const double b01 = Nu[0][i] * Nv[1][j];
			const double b20 = Nu[2][i] * Nv[0][j];
			const double b11 = Nu[1][i] * Nv[1][j];
			const double b02 = Nu[0][i] * Nv[2][j];

			A += Pw * b00;
			Au += Pw * b10;
			Av += Pw * b01;
			Auu += Pw * b20;
			Auv += Pw * b11;
			Avv += Pw * b02;

			W += ww * b00;
			Wu += ww * b10;
			Wv += ww * b01;
			Wuu += ww * b20;
			Wuv += ww * b11;
			Wvv += ww * b02;
		}
	}

	if (std::fabs(W) < 1.e-14)
		return;

	const double W2 = W * W;
	const double W3 = W2 * W;

	du = (Au * W - A * Wu) / W2;
	dv = (Av * W - A * Wv) / W2;

	duu = (Auu * W2 - A * (W * Wuu) - Au * (2. * W * Wu) + A * (2. * Wu * Wu)) / W3;
	dvv = (Avv * W2 - A * (W * Wvv) - Av * (2. * W * Wv) + A * (2. * Wv * Wv)) / W3;
	duv = (Auv * W2 - Au * (W * Wv) - Av * (W * Wu) - A * (W * Wuv) + A * (2. * Wu * Wv)) / W3;
}

bool NurbsSurface::normal(double u, double v, Point3& n) const
{
	Point3 du, dv, duu, duv, dvv;
	evaluate_derivatives(u, v, du, dv, duu, duv, dvv);

	n = du.cross_product(dv);
	const double norm2 = n.norm_square();
	if (norm2 < 1.e-20)
	{
		n = Point3();
		return false;
	}
	n /= std::sqrt(norm2);
	return true;
}

bool NurbsSurface::curvature(double u, double v, double& gaussian, double& mean) const
{
	gaussian = 0.;
	mean = 0.;

	Point3 du, dv, duu, duv, dvv;
	evaluate_derivatives(u, v, du, dv, duu, duv, dvv);

	Point3 n = du.cross_product(dv);
	const double n2 = n.norm_square();
	if (n2 < 1.e-20)
		return false;
	n /= std::sqrt(n2);

	const double E = du.dot_product(du);
	const double F = du.dot_product(dv);
	const double G = dv.dot_product(dv);

	const double L = n.dot_product(duu);
	const double M = n.dot_product(duv);
	const double N = n.dot_product(dvv);

	const double den = E * G - F * F;
	if (std::fabs(den) < 1.e-20)
		return false;

	gaussian = (L * N - M * M) / den;
	mean = (E * N - 2. * F * M + G * L) / (2. * den);
	return true;
}

void NurbsSurface::evaluate_clamped(double u, double v, Point3& p) const
{
	evaluate(std::clamp(u, 0.0, 1.0), std::clamp(v, 0.0, 1.0), p);
}

void NurbsSurface::evaluate_partials(double u, double v, Point3& du, Point3& dv) const
{
	const double h = 1.e-4;

	double u1 = std::clamp(u - h, 0.0, 1.0);
	double u2 = std::clamp(u + h, 0.0, 1.0);
	double v1 = std::clamp(v - h, 0.0, 1.0);
	double v2 = std::clamp(v + h, 0.0, 1.0);

	Point3 pu1; evaluate_clamped(u1, v, pu1);
	Point3 pu2; evaluate_clamped(u2, v, pu2);
	Point3 pv1; evaluate_clamped(u, v1, pv1);
	Point3 pv2; evaluate_clamped(u, v2, pv2);

	double duStep = (u2 - u1);
	double dvStep = (v2 - v1);

	if (duStep <= 0.)
		du = Point3();
	else
		du = (pu2 - pu1) / duStep;

	if (dvStep <= 0.)
		dv = Point3();
	else
		dv = (pv2 - pv1) / dvStep;
}



void NurbsSurface::project_point_on_surface(const Point3& target, double& u, double& v, Point3& projected) const
{
	u = std::clamp(u, 0.0, 1.0);
	v = std::clamp(v, 0.0, 1.0);

	for (int i = 0; i < 10; ++i)
	{
		Point3 p;
		evaluate(u, v, p);
		Point3 r = p - target;

		Point3 du, dv;
		evaluate_partials(u, v, du, dv);

		double a11 = du.dot_product(du) + 1.e-12;
		double a12 = du.dot_product(dv);
		double a22 = dv.dot_product(dv) + 1.e-12;
		double b1 = du.dot_product(r);
		double b2 = dv.dot_product(r);

		double det = a11 * a22 - a12 * a12;
		if (std::fabs(det) < 1.e-18)
			break;

		double dU = (-b1 * a22 + b2 * a12) / det;
		double dV = (-a11 * b2 + a12 * b1) / det;

		u = std::clamp(u + dU, 0.0, 1.0);
		v = std::clamp(v + dV, 0.0, 1.0);

		if (dU * dU + dV * dV < 1.e-14)
			break;
	}

	evaluate(u, v, projected);
}

bool NurbsSurface::is_trimmed() const
{
	return false;
}

void NurbsSurface::reverse_u()
{
	const int nU = nb_points_u();
	const int nV = nb_points_v();
	if (nU <= 0 || nV <= 0)
		return;

	std::vector<Point3> newPoints(points().size());
	std::vector<double> newWeights(weights().size());

	for (int v = 0; v < nV; ++v)
		for (int u = 0; u < nU; ++u)
		{
			const int srcIdx = u + nU * v;
			const int dstIdx = (nU - 1 - u) + nU * v;
			newPoints[dstIdx] = points()[srcIdx];
			newWeights[dstIdx] = weights()[srcIdx];
		}

	std::vector<double> newKnotsU = knots_u();
	if (!newKnotsU.empty())
	{
		const int m = (int)newKnotsU.size();
		for (int i = 0; i < m; ++i)
			newKnotsU[i] = 1. - knots_u()[m - 1 - i];
	}

	set_points(newPoints, nU, nV);
	set_weights(newWeights);
	set_knots_u(newKnotsU);
}

NurbsSurface NurbsSurface::reversed_u() const
{
	NurbsSurface copy = *this;
	copy.reverse_u();
	return copy;
}
