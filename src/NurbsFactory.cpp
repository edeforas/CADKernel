#include "NurbsFactory.h"

#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "NurbsSolid.h"
#include "NurbsRuled.h"
#include "NurbsRevolve.h"
#include "NurbsUtil.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_circle(double dRadius, NurbsCurve& nc)
{
	nc.clear();

	NurbsCurve n;
	int degree = 2;
	std::vector<double> knots = { 0., 0., 0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1., 1., 1. };
	std::vector<double> weights = { 1.,1. / std::sqrt(2.),1.,1. / std::sqrt(2.),1.,1. / std::sqrt(2.),1.,1. / std::sqrt(2.),1. };
	std::vector<Point3> points = {
		Point3(1.,0.,0.) * dRadius,Point3(1.,1.,0.) * dRadius,Point3(0.,1.,0.) * dRadius,
		Point3(-1.,1.,0.) * dRadius,Point3(-1.,0.,0.) * dRadius,Point3(-1.,-1.,0.) * dRadius,
		Point3(0.,-1.,0.) * dRadius,Point3(1.,-1.,0.) * dRadius,Point3(1.,0.,0.) * dRadius
	};

	nc.set_degree(degree);
	nc.set_knots(knots);
	nc.set_weights(weights);
	nc.set_points(points);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_disk(double dRadius, NurbsSurface& ns)
{
	ns.clear();

	//create profile curve
	NurbsCurve nc;
	std::vector<Point3> points = {
		Point3(0., 0.,0.),
		Point3(dRadius,0.,0.),
	};

	nc.set_degree(1);
	nc.set_points(points);
	nc.set_uniform();
	nc.set_equals_weights();

	NurbsRevolve nr;
	nr.revolve(nc, ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_filled_surface(const NurbsCurve& nc, NurbsSurface& ns)
{
	ns.clear();

	if (nc.nb_points() < 2)
		return;

	if (nc.weights().size() != nc.nb_points())
		return;

	const std::vector<Point3>& curvePoints = nc.points();
	const std::vector<double>& curveWeights = nc.weights();

	int nbCenterSamples = nc.nb_points();
	if (nc.is_closed() && nc.nb_points() > 2)
		nbCenterSamples = nc.nb_points() - 1;

	Point3 center;
	double totalWeight = 0.;
	for (int i = 0; i < nbCenterSamples; i++)
	{
		double w = curveWeights[i];
		if (std::fabs(w) < 1.e-12)
			continue;

		center += curvePoints[i] * w;
		totalWeight += w;
	}

	if (std::fabs(totalWeight) < 1.e-12)
		center = curvePoints[0];
	else
		center /= totalWeight;

	std::vector<Point3> surfacePoints = curvePoints;
	std::vector<double> surfaceWeights = curveWeights;

	for (int i = 0; i < nc.nb_points(); i++)
	{
		surfacePoints.push_back(center);
		surfaceWeights.push_back(curveWeights[i]);
	}

	ns.set_degree(nc.degree(), 1);
	ns.set_knots_u(nc.knots());
	ns.set_knots_v({ 0., 0., 1., 1. });
	ns.set_points(surfacePoints, nc.nb_points(), 2);
	ns.set_weights(surfaceWeights);
	ns.set_closed_u(nc.is_closed());
	ns.set_closed_v(false);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_sphere(double dRadius,NurbsSurface& ns)
{
	ns.clear();

	//create profile curve: half circle
	NurbsCurve nc;
	std::vector<Point3> points = {
		Point3(0., 0.,-1.)* dRadius,
		Point3(0.,1.,-1.)* dRadius,
		Point3(0., 1.,0.)* dRadius,
		Point3(0.,1.,1.)* dRadius,
		Point3(0.,0.,1.)* dRadius
	};

	std::vector<double> weights = { 1. , 1. / std::sqrt(2.),1., 1. / std::sqrt(2.), 1. };
	std::vector<double> knots = { 0., 0., 0., 0.5, 0.5, 1., 1., 1. };
	nc.set_degree(2);
	nc.set_points(points);
	nc.set_weights(weights);
	nc.set_knots(knots);

	NurbsRevolve nr;
	nr.revolve(nc, ns);

	for (auto& p : ns.points())
	{
		p.y() = -p.y();
	}

}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_torus(double dMajorRadius,double dMinorRadius, NurbsSurface& ns)
{
	ns.clear();

	//create profile curve: circle and transform it
	NurbsCurve nc;
	NurbsFactory::create_circle(dMinorRadius, nc);
	for (int i = 0; i < nc.nb_points(); i++)
	{
		Point3& p = nc.points()[i];

		p.x() += dMajorRadius;

		double dz = p.z();
		p.z() = p.y();
		p.y() = -dz;
	}

	// revolve
	NurbsRevolve nr;
	nr.revolve(nc, ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_cylinder(double dRadius, double dHeight, NurbsSolid& n)
{
	n.clear();

	NurbsSurface s1, s2, s3;

	create_disk(dRadius, s1);
	for (auto& p : s1.points())
	{
		p.z() -= dHeight / 2.;
		p.y() = -p.y();
	}

	create_disk(dRadius, s3);
	for (auto& p : s3.points())
	{
		p.z() += dHeight / 2.;
	}

	//create profile curve
	NurbsCurve nc;
	std::vector<Point3> points = {
		Point3(dRadius,0.,-dHeight / 2.),Point3(dRadius,0.,dHeight / 2.)
	};

	NurbsUtil::create_curve_from_points(points, 1, nc);
	NurbsRevolve nr;
	nr.revolve(nc, s2);

	for (auto& p : s2.points())
	{
		p.y() = -p.y();
	}

	n.add_surface(s1);
	n.add_surface(s2);
	n.add_surface(s3);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_torus(double dMajorRadius, double dMinorRadius, NurbsSolid& n)
{
	n.clear();
	NurbsSurface ns;
	create_torus(dMajorRadius, dMinorRadius, ns);

	for (auto& p : ns.points())
	{
		p.y() = -p.y();
	}

	n.add_surface(ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_sphere(double dRadius, NurbsSolid& n)
{
	n.clear();
	NurbsSurface ns;
	create_sphere(dRadius, ns);
	n.add_surface(ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_triangle(const Point3& p1, const Point3& p2, const Point3& p3, NurbsSurface& ns)
{
	create_quad(p1, p2, p3, p3, ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_quad(const Point3& p1, const Point3& p2, const Point3& p3, const Point3& p4, NurbsSurface& ns)
{
	NurbsCurve n1,n2;
	
	NurbsUtil::create_curve_from_points({ p1,p2 }, 1, n1);
	NurbsUtil::create_curve_from_points({ p4,p3 }, 1, n2);
	NurbsRuled nr;
	nr.create_ruled(n1, n2, ns);
}
///////////////////////////////////////////////////////////////////////////
void NurbsFactory::create_box(double dSizeX, double dSizeY, double dSizeZ, NurbsSolid& n)
{
	n.clear();

	NurbsSurface ns;

	// +Z face
	Point3 p1(dSizeX / 2., dSizeY / 2., dSizeZ / 2.);
	Point3 p2(dSizeX / 2., -dSizeY / 2., dSizeZ / 2.);
	Point3 p3(-dSizeX / 2., -dSizeY / 2., dSizeZ / 2.);
	Point3 p4(-dSizeX / 2., dSizeY / 2., dSizeZ / 2.);

	// -Z face
	Point3 p5(dSizeX / 2., dSizeY / 2., -dSizeZ / 2.);
	Point3 p6(dSizeX / 2., -dSizeY / 2., -dSizeZ / 2.);
	Point3 p7(-dSizeX / 2., -dSizeY / 2., -dSizeZ / 2.);
	Point3 p8(-dSizeX / 2., dSizeY / 2., -dSizeZ / 2.);

	create_quad(p1, p4, p3, p2, ns);
	n.add_surface(ns);

	create_quad(p5, p6, p7, p8, ns);
	n.add_surface(ns);

	create_quad(p1, p2, p6, p5, ns);
	n.add_surface(ns);

	create_quad(p3, p4, p8, p7, ns);
	n.add_surface(ns);

	create_quad(p3, p7, p6, p2, ns);
	n.add_surface(ns);

	create_quad(p1, p5, p8, p4, ns);
	n.add_surface(ns);
}
///////////////////////////////////////////////////////////////////////////


