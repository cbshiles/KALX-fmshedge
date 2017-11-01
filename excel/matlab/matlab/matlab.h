// matlab.h - header file for project
// Uncomment the following line to use features for Excel2007 and above.
//#define EXCEL12
#include "ooura.h"
#include "normal.h"
#include "xll/xll.h"
#include "math.h"
#include <random>
#include <functional>
#include "secant.hpp"
#include <tuple>

#ifndef CATEGORY
#define CATEGORY _T("MATLAB")
#endif

typedef xll::traits<XLOPERX>::xcstr xcstr; // pointer to const string
typedef xll::traits<XLOPERX>::xword xword; // use for OPER and FP indices
typedef xll::traits<XLOPERX>::xfp xfp; // use for OPER and FP indices

inline
double normcdf(double z)
{
	return normal_cdf(z);
}

/*Calculate the price of a call given parameters
@S: Initial Price of underlying
@K: Strike of the Call
@r: Risk-Free Interest Rate
@sigma: historical volatility of underlying asset
@tMaturity: tMaturity to Maturity of the call
*/
inline double blsprice(double S,double K,double r,double sigma,double t){

double t_sqrt = sqrt(t);
double d1 = (log(S/K)+r*t)/(sigma*t_sqrt)+0.5*sigma*t_sqrt;
double d2 = d1-(sigma*t_sqrt);
double price  = S*normal_cdf(d1)-K*exp(-r*t)*normal_cdf(d2);

return price;

};

/*Calculate the delta of a call given parameters
@S: Initial Price of underlying
@K: Strike of the Call
@r: Risk-Free Interest Rate
@sigma: historical volatility of underlying asset
@t: t to Maturity of the call
*/
inline double blsdelta(double  S,double  K,double  r,double  sigma,double  t){

double t_sqrt = sqrt(t);
double d1 = (log(S/K)+r*t)/(sigma*t_sqrt)+0.5*sigma*t_sqrt;
double delta = normal_cdf(d1);

return delta;
};

std::tuple<double,bool> hitting_times_bisection(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA);
std::tuple<double,double> median_hitting_time(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA);