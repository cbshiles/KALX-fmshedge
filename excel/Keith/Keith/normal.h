// normal.h - Normally distributed random variables.
#pragma once
#include <cmath>
#include <limits>
#include "ooura.h"

double normal_pdf(double);

template <class T>
double normal_cdf(double);

template <class T>
double normal_inv(double);

#ifndef M_PI
#define M_PI		3.1415926535897931e+00
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT2PI
#define M_SQRT2PI sqrt(2*M_PI)
#endif

inline
double normal_pdf(double x)
{
	return exp(-x*x/2)/M_SQRT2PI;
}

template<class T> inline double normal_cdf(double);

// From: bill.daly@tradition-ny.com (Bill Daly)
// Message-ID: <4fa0d5e.0402241117.53c8da2a@posting.google.com>
// Message-ID: <4fa0d5e.0402181449.6dd203db@posting.google.com>

class daly {};
template <> inline double
normal_cdf<daly>(double x)
{
	long double s=x,t=0,b=x,x2=x*x,i=1;
	
	while(s!=t)
		s=(t=s)+(b*=x2/(i+=2));

	return  static_cast<double>(.5+s*exp(-.5*x2-.91893853320467274178L));
}

template<> inline double
normal_cdf<ooura>(double z)
{
	return 1 - derfc<ooura>(z/M_SQRT2)/2;
}

template<> inline double
normal_inv<ooura>(double p)
{
	return M_SQRT2*dierfc<ooura>(2*(1 - p));
}

// Abramowitz and Stegun. Handbook of Mathematical functions 26.24
class AS_P1 {};
template<> inline double
normal_cdf<AS_P1>(double x)
{
	return x == 0 ? 0.5
		: x > 0 ? 0.5 * (1 + sqrt(1 - exp(-2*x*x/M_PI)))
		: 0.5 * (1 - sqrt(1 - exp(-2*x*x/M_PI)));
}
class AS_P2 {};
template<> inline double
normal_cdf<AS_P2>(double x)
{
	return x > 1.4 ? 1 - normal_pdf(x)*(sqrt(4 + x*x) - x)/2
		: x < -1.4 ? normal_pdf(x)*(sqrt(4 + x*x) + x)/2
		: normal_cdf<AS_P1>(x);
}

template<> inline double
normal_inv<AS_P1>(double p)
{
	return p > 0.5
		?  sqrt(-(M_PI/2)*log(4*p*(1-p)))
		: -sqrt(-(M_PI/2)*log(4*p*(1-p)));
}

class zyang {};
template<> inline double
normal_inv<zyang>(double p)
{
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0;

   if (p > 0.5)
   	return -normal_inv<zyang>(1 - p);

   if (p < 1e-20)
   	return std::numeric_limits<double>::quiet_NaN();

   y = sqrt (log(1/(p*p)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);

   return  -z;
}

template<> inline double
normal_inv<daly>(double p)
{
	double x, y = 0;

	if (p < 0.5)
		return -normal_inv<daly>(1 - p);

	for (x = normal_inv<AS_P1>(p); x != y; x = y) {
		y = x - (normal_cdf<daly>(x) - p)/normal_pdf(x); 
	}

	return x;
}
