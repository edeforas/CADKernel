#ifndef _Geometry_
#define _Geometry_

class Plane3;

///////////////////////////////////////////////////////////////////////////
class Point3
{
public:
	Point3();
    Point3(double x, double y, double z);

    Point3& operator+=(const Point3& p);
	Point3& operator-=(const Point3& p);
	Point3& operator*=(double d);
	Point3& operator/=(double d);
	Point3 operator+(const Point3& p) const;
	Point3 operator*(double d) const;
	Point3 operator/(double d) const;
    Point3 operator-(const Point3& p) const;

	double& x();
	double& y();
	double& z();

	const double& x() const;
	const double& y() const;
	const double& z() const;

    double distance_square(const Point3& p) const;
    double dot_product(const Point3& p) const;
    Point3 cross_product(const Point3& p) const;
	double norm() const;
	double norm_square() const;
	void normalize();
	Point3 normalized() const;
	void sanitize();

    double _x, _y, _z;
};
///////////////////////////////////////////////////////////////////////////
class Line3
{
public:
	Line3();
	Line3(const Point3& p1, const Point3& p2);

	const Point3& p1() const;
	const Point3& p2() const;

	void set_p1(const Point3& p);
	void set_p2(const Point3& p);

private:
	Point3 _p1;
	Point3 _p2;
};
///////////////////////////////////////////////////////////////////////////
class Segment3
{
public:
	Segment3();
	Segment3(const Point3& p1, const Point3& p2);

	const Point3& p1() const;
	const Point3& p2() const;

	void set_p1(const Point3& p);
	void set_p2(const Point3& p);

	double norm() const;
	double norm_square() const;

	bool intersect(const Segment3& s, Point3& pIntersection) const;

private:
	Point3 _p1;
	Point3 _p2;
};
///////////////////////////////////////////////////////////////////////////
class Triangle3
{
public:
    Triangle3();
    Triangle3(const Point3& p1, const Point3& p2, const Point3& p3);
	
	const Point3& p1() const;
	const Point3& p2() const;
	const Point3& p3() const;

	void set_p1(const Point3& p);
	void set_p2(const Point3& p);
	void set_p3(const Point3& p);

	bool cutted_by(const Plane3& p) const;
	bool contains(const Point3& p) const;
	bool intersect_with(const Segment3& s, Point3& pIntersection) const;
	bool intersect_with(const Triangle3& t) const;

	double surface() const;
	Point3 normal() const; // normal vector (normalized)
	Point3 orthogonal() const; // compute orthogonal vector (normal to triangle but not normalized())

private:
	Point3 _p1, _p2, _p3;
};
///////////////////////////////////////////////////////////////////////////
class Plane3
{
public:
	Plane3();
	Plane3(const Point3& p1, const Point3& p2, const Point3& p3);
	explicit Plane3(const Triangle3& t);

	void compute_using(const Point3& p1, const Point3& p2, const Point3& p3);
	void compute_using(const Triangle3& t);
	
	double distance_to(const Point3& p) const;

	void get(double& a, double& b, double& c, double& d) const;
	
	bool intersect_with(const Segment3& s, Point3& pIntersection) const;
	bool intersect_with(const Plane3 &p, const Line3& pIntersection) const;

	void project_point(const Point3& p, Point3& projected) const;

	Point3 normal() const; // normal vector (normalized)

private:
		double _a, _b, _c, _d;
};
///////////////////////////////////////////////////////////////////////////
class BoundingBox3 // AABB: Axis Aligned Bounding Box
{
public:
	BoundingBox3();
	explicit BoundingBox3(const Triangle3 & t);

	void add(const Point3& p1);
	bool intersect_with(const BoundingBox3& b) const;

private:
	double _x1, _y1, _z1;
	double _x2, _y2, _z2;
	bool _bInitialized;
};

///////////////////////////////////////////////////////////////////////////
#endif