#include "BezierUtil.h"

double BezierUtil::square(double a)
{
	return a * a;
}

double BezierUtil::cubic(double a)
{
	return a * a * a;
}

double BezierUtil::Bernstein01(double u)
{
	return 1. - u;
}

double BezierUtil::Bernstein11(double u)
{
	return u;
}

double BezierUtil::Bernstein02(double u)
{
	return square(1. - u);
}

double BezierUtil::Bernstein12(double u)
{
	return 2. * u * (1. - u);
}

double BezierUtil::Bernstein22(double u)
{
	return square(u);
}

double BezierUtil::Bernstein03(double u)
{
	return cubic(1. - u);
}

double BezierUtil::Bernstein13(double u)
{
	return 3. * u * square(1. - u);
}

double BezierUtil::Bernstein23(double u)
{
	return 3. * square(u) * (1. - u);
}

double BezierUtil::Bernstein33(double u)
{
	return cubic(u);
}
