// header.h - header file for project
// Uncomment the following line to use features for Excel2007 and above.
//#define EXCEL12
#pragma once
#include "xll/xll.h"
#include "normal.h"
#include "ooura.h"

//#include "hedge.h"



#ifndef CATEGORY
//Function Wizard Category
#define CATEGORY _T("Cornell_MEng_Project_KALX")
#endif

typedef xll::traits<XLOPERX>::xcstr xcstr; // pointer to const string
typedef xll::traits<XLOPERX>::xword xword; // use for OPER and FP indices

//inline: Define once 

/*Calculate the price of a call given parameters
@S: Initial Price of underlying
@K: Strike of the Call
@r: Risk-Free Interest Rate
@sigma: historical volatility of underlying asset
@tMaturity: tMaturity to Maturity of the call
*/
inline double price_call(double S,
	     double K,
	     double r,
	     double sigma,
	     double t
	     ){

//Black-Scholes Formula
double t_sqrt = sqrt(t);
double d1 = (log(S/K)+r*t)/(sigma*t_sqrt)+0.5*sigma*t_sqrt;
double d2 = d1-(sigma*t_sqrt);
double price  = S*normal_cdf(d1)-K*exp(-r*t)*normal_cdf(d2);

return price;
}


/*Calculate the delta of a call given parameters
@S: Initial Price of underlying
@K: Strike of the Call
@r: Risk-Free Interest Rate
@sigma: historical volatility of underlying asset
@t: t to Maturity of the call
*/
inline double delta_call(double  S,
	     double  K,
	     double  r,
	     double  sigma,
	     double  t){

double t_sqrt = sqrt(t);
double d1 = (log(S/K)+r*t)/(sigma*t_sqrt)+0.5*sigma*t_sqrt;

double delta = normal_cdf(d1);

return delta;

}

