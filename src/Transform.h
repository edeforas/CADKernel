#ifndef _Transform_
#define _Transform_

#include "Geometry.h"
#include <vector>
using namespace std;

class Transform
{
public:
	Transform();

	virtual void apply(Point3& p) const = 0;
	void apply_all(vector<Point3>& vp) const;
};
///////////////////////////////////////////////////////////////////////////
class Translation : public Transform
{
public:
	Translation(const Point3& translation);
	Translation(double tx, double ty, double tz);

	virtual void apply(Point3& p) const override;

private:
	Point3 _translation;
};
///////////////////////////////////////////////////////////////////////////
class Scale : public Transform
{
public:
	Scale(double scaleX = 1., double scaleY = 1., double scaleZ = 1.);

	virtual void apply(Point3& p) const override;

private:
	double 	_scaleX, _scaleY, _scaleZ;
};
///////////////////////////////////////////////////////////////////////////
class RotationEulerAngles : public Transform
{
public:
	RotationEulerAngles(double yaw = 0., double pitch = 0., double roll = 0.);
	void set_angles(double yaw = 0., double pitch = 0., double roll = 0.);
	virtual void apply(Point3& p) const override;
private:
	double 	_yaw, _pitch, _roll;
	double _yawCos, _yawSin, _pitchCos, _pitchSin, _rollCos, _rollSin;
};
///////////////////////////////////////////////////////////////////////////
class RotationQuaternion : public Transform
{
public:
	RotationQuaternion();
	RotationQuaternion(double q[4]);
	void set_quaternion(double q[4]);
	void from_angle_axis(double dAngle, const Point3& axis);
	virtual void apply(Point3& p) const override;
private:
	double _q[4];
};
///////////////////////////////////////////////////////////////////////////
class RotationMatrix : public Transform
{
public:
	RotationMatrix();
	RotationMatrix(double m[9]);
	void set_matrix(double m[9]);
	void get_matrix(double m[9]) const;
	void from_angle_axis(double dAngle, const Point3& axis);
	virtual void apply(Point3& p) const override;
private:
	double _m[9];
};
///////////////////////////////////////////////////////////////////////////
class RotationAngleAxis : public Transform
{
public:
	RotationAngleAxis();
	RotationAngleAxis(double dAngle, const Point3& axis);
	void set_angle_axis(double dAngle, const Point3& axis);
	void from_angle_axis(double dAngle, const Point3& axis);
	virtual void apply(Point3& p) const override;
private:
	RotationMatrix _rot;
};
///////////////////////////////////////////////////////////////////////////
#endif