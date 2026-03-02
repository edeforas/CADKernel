#include "Transform.h"

#include <cmath>
#include <cstring> //for memcpy

///////////////////////////////////////////////////////////////////////////
Transform::Transform()
{
}
void Transform::apply_all(vector<Point3>& vp) const
{
	for (auto& p : vp)
		apply(p);
}
///////////////////////////////////////////////////////////////////////////
Translation::Translation(const Point3& translation):
	_translation(translation)
{ }

Translation::Translation(double tx, double ty, double tz):
	_translation(Point3(tx, ty, tz))
{ }

void Translation::apply(Point3& p) const
{
	p += _translation;
}
///////////////////////////////////////////////////////////////////////////
Scale::Scale(double scaleX, double scaleY, double scaleZ)
{
	_scaleX = scaleX;
	_scaleY = scaleY;
	_scaleZ = scaleZ;
}
void Scale::apply(Point3& p) const
{
	p.x() *= _scaleX;
	p.y() *= _scaleY;
	p.z() *= _scaleZ;
}
///////////////////////////////////////////////////////////////////////////
RotationEulerAngles::RotationEulerAngles(double yaw, double pitch, double roll)
{
	set_angles(yaw, pitch, roll);
}
void RotationEulerAngles::set_angles(double yaw, double pitch, double roll)
{
	double dDegToRad = 2. * 3.14159265359 / 360.;

	_yaw = yaw;
	_pitch = pitch;
	_roll = roll;

	_yawCos = cos(_yaw * dDegToRad);
	_yawSin = sin(_yaw * dDegToRad);
	_pitchCos = cos(_pitch * dDegToRad);
	_pitchSin = sin(_pitch * dDegToRad);
	_rollCos = cos(_roll * dDegToRad);
	_rollSin = sin(_roll * dDegToRad);
}
void RotationEulerAngles::apply(Point3& p) const
{
	//Rotate Yaw then Pitch , then Roll

	//yaw rotation
	double tmp = p.x();
	p.x() = p.x() * _yawCos + p.z() * _yawSin;
	p.z() = p.z() * _yawCos - tmp * _yawSin;

	//pitch rotation
	tmp = p.y();
	p.y() = p.z() * _pitchSin - p.y() * _pitchCos;
	p.z() = p.z() * _pitchCos + tmp * _pitchSin;

	//roll rotation
	tmp = p.x();
	p.x() = p.x() * _rollCos + p.y() * _rollSin;
	p.y() = p.y() * _rollCos - tmp * _rollSin;
}
///////////////////////////////////////////////////////////////////////////
RotationQuaternion::RotationQuaternion()
{
	_q[0] = 1.;
	_q[1] = 0.;
	_q[2] = 0.;
	_q[3] = 0.;
}
RotationQuaternion::RotationQuaternion(double q[4])
{
	set_quaternion(q);
}
void RotationQuaternion::set_quaternion(double q[4])
{
	std::memcpy(_q, q, 4 * sizeof(double));
}
void RotationQuaternion::from_angle_axis(double dAngle, const Point3& axis)
{
	double dDegToRad = 2. * 3.14159265359 / 360.;
	double s = sin(dAngle * dDegToRad * 0.5);
	double c = cos(dAngle * dDegToRad * 0.5);

	_q[0] = c;
	_q[1] = axis.x() * s;
	_q[2] = axis.y() * s;
	_q[3] = axis.z() * s;
}

void RotationQuaternion::apply(Point3& p) const
{
	// from https://fr.wikipedia.org/wiki/Quaternions_et_rotation_dans_l%27espace

	double a = _q[0];
	double b = _q[1];
	double c = _q[2];
	double d = _q[3];

	double t2 = a * b;
	double t3 = a * c;
	double t4 = a * d;
	double t5 = -b * b;
	double t6 = b * c;
	double t7 = b * d;
	double t8 = -c * c;
	double t9 = c * d;
	double t10 = -d * d;

	double v1 = p.x();
	double v2 = p.y();
	double v3 = p.z();

	p.x() = 2. * ((t8 + t10) * v1 + (t6 - t4) * v2 + (t3 + t7) * v3) + v1;
	p.y() = 2. * ((t4 + t6) * v1 + (t5 + t10) * v2 + (t9 - t2) * v3) + v2;
	p.z() = 2. * ((t7 - t3) * v1 + (t2 + t9) * v2 + (t5 + t8) * v3) + v3;
}
///////////////////////////////////////////////////////////////////////////
RotationMatrix::RotationMatrix()
{
	_m[0] = 1.;
	_m[1] = 0.;
	_m[2] = 0.;
	_m[3] = 0.;
	_m[4] = 1.;
	_m[5] = 0.;
	_m[6] = 0.;
	_m[7] = 0.;
	_m[8] = 1.;
}
RotationMatrix::RotationMatrix(double m[9])
{
	set_matrix(m);
}
void RotationMatrix::set_matrix(double m[9])
{
	std::memcpy(_m, m, 9 * sizeof(double));
}

void RotationMatrix::get_matrix(double m[9]) const
{
	std::memcpy(m, _m, 9 * sizeof(double));
}

void RotationMatrix::from_angle_axis(double dAngle, const Point3& axis)
{
	// from https://fr.wikipedia.org/wiki/Quaternions_et_rotation_dans_l%27espace

	double v1 = axis.x();
	double v2 = axis.y();
	double v3 = axis.z();

	double dDegToRad = 2. * 3.14159265359 / 360.;
	double t1 = cos(dAngle * dDegToRad);
	double t2 = 1. - t1;
	double t3 = v1 * v1;
	double t6 = t2 * v1;
	double t7 = t6 * v2;
	double t8 = sin(dAngle * dDegToRad);
	double t9 = t8 * v3;
	double t11 = t6 * v3;
	double t12 = t8 * v2;
	double t15 = v2 * v2;
	double t19 = t2 * v2 * v3;
	double t20 = t8 * v1;
	double t24 = v3 * v3;

	_m[0] = t1 + t2 * t3;
	_m[1] = t7 - t9;
	_m[2] = t11 + t12;
	_m[3] = t7 + t9;
	_m[4] = t1 + t2 * t15;
	_m[5] = t19 - t20;
	_m[6] = t11 - t12;
	_m[7] = t19 + t20;
	_m[8] = t1 + t2 * t24;
}
void RotationMatrix::apply(Point3& p) const
{
	//left matrix x'=m*x -> p is a column vector
	double x = p.x();
	double y = p.y();
	double z = p.z();

	p.x() = x * _m[0] + y * _m[1] + z * _m[2];
	p.y() = x * _m[3] + y * _m[4] + z * _m[5];
	p.z() = x * _m[6] + y * _m[7] + z * _m[8];
}
///////////////////////////////////////////////////////////////////////////
RotationAngleAxis::RotationAngleAxis()
{
}
RotationAngleAxis::RotationAngleAxis(double dAngle, const Point3& axis)
{
	_rot.from_angle_axis(dAngle, axis);
}
void RotationAngleAxis::set_angle_axis(double dAngle, const Point3& axis)
{
	_rot.from_angle_axis(dAngle, axis);
}
void RotationAngleAxis::from_angle_axis(double dAngle, const Point3& axis)
{
	_rot.from_angle_axis(dAngle, axis);
}
void RotationAngleAxis::apply(Point3& p) const
{
	_rot.apply(p);
}
///////////////////////////////////////////////////////////////////////////
