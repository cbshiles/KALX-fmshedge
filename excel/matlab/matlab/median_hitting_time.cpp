#include "matlab.h"

//auto rand = std::bind( std::uniform_real_distribution<double>(0,1), std::default_random_engine());

using namespace xll;

//Function to compute random hitting times for a geometric BM with two barrier levels (upper and lower)
//This function returns to a 1 by 2 matrix tuple with upper median hitting times for the first column and the lower median hitting times for the second column 

//Inputs
//S_CURRENT-Current value of the spot price (scalar)
//S_UPPER-Current level of the upper barrier (scalar)
//S_LOWER-Current level of the lower barrier (scalar)
//R - Risk free interest rate (drift for the risk neutral probability measure)
//SIGMA - Volatility (difusion term)

//Outputs
//MEDIAN_HIT_TIME_UPPER - median of the hitting time at the upper barrier
//MEDIAN_HIT_TIME_LOWER - median of the hitting time at the lower barrier

std::tuple<double,double> median_hitting_time(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA)
{
//Algorithm
double year2hour = 12*21*24; // Rescale the time dimension to hours
R = R/year2hour; // rescale R to hours
SIGMA = SIGMA/sqrt(year2hour); // Rescale SIGMA to hours

double alpha = (1/SIGMA)*log(S_UPPER/S_CURRENT); // Defines the equivalent level of the upper barrier of the GBM to the standard BM
double beta = (1/SIGMA)*log(S_LOWER/S_CURRENT); // Defines the equivalent level of the upper barrier of the GBM to the standard BM
double nu = (R/SIGMA) - (SIGMA/2); // Defines parameter nu

//Generates 0.5 and 0.5 as the uniform "random numbers" to ensure we can compute the median of both.
double U1 = .5; 
double U2 = .5; 

//Calls the bisection method for both upper and lower barriers
double MEDIAN_HIT_TIME_UPPER = root1d::secant_bisect([=](double t){return U1 - normcdf((-alpha + nu*t)/sqrt(t)) - exp(2*nu*alpha)*normcdf((-nu*t - alpha)/sqrt(t));},400.,800.,1e-11,100,0.,10000.);
double MEDIAN_HIT_TIME_LOWER = root1d::secant_bisect([=](double t){return U2 - normcdf((beta - nu*t)/sqrt(t)) - exp(2*nu*beta)*normcdf((nu*t + beta)/sqrt(t));},400.,800.,1e-11,100,0.,10000.);

return std::make_tuple(MEDIAN_HIT_TIME_UPPER,MEDIAN_HIT_TIME_LOWER);
}

static AddInX xai_median_hitting_time(
	// Return type, C name, Excel name.
	FunctionX(XLL_LPOPERX, _T("?xll_median_hitting_time"), _T("XLL.median_hitting_time"))
	.Arg(XLL_DOUBLEX, _T("S_CURRENT"), _T("Current stock price. (scalar) "),100)
	.Arg(XLL_DOUBLEX, _T("S_UPPER"), _T("Upper bound of limit. (scalar)"),105)
	.Arg(XLL_DOUBLEX, _T("S_LOWER"), _T("is the level of current price in stocks. (scalar) "),95)
	.Arg(XLL_DOUBLEX, _T("R"), _T("is the normal R.V w/ Mean 0 and Std 1. (scalar) "),0.01)
	.Arg(XLL_DOUBLEX, _T("SIGMA"), _T("is the level of volatility of the time step (dt). (scalar) "),.25)
	.Category(CATEGORY)
	.FunctionHelp(_T("Generates the change in stock price, according to BlackScholes PDE."))
	.Documentation(_T("Optional documentation for function: None."))
);
LPOPERX WINAPI xll_median_hitting_time(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA)
{
	#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	static OPERX o(2,1);

	// function body goes here
	try{
	std::tuple<double,double> db = median_hitting_time(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA);
	o[0] = std::get<0>(db);
	o[1] = std::get<1>(db);	
	}
	catch(
	const std::exception& ex
	)
	{
		XLL_ERROR (ex.what());
	}

	return &o;
}
