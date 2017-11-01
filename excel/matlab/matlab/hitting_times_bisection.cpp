#include "matlab.h"



auto rand1 = std::bind( std::uniform_real_distribution<double>(0,1), std::default_random_engine());



using namespace xll;



//Function to compute random hitting times for a GEOMETRIC BM with two barrrier levels (upper and lower)

//This function returns a 1 by 2 tuple with the hitting time for the first column and the flag variable for the second column.



//Inputs

//S_CURRENT - Current value of the spot price (scalar)

//S_UPPER - Current level of the upper barrier (scalar)

//S_LOWER - Current level of the lower barrier (scalar)

//R - Risk free interest rate (drift for the risk neutral probability measure)

//SIGMA - Volatility (difusion term)



//Outputs

//HIT_TIME - hitting time for either upper or lower levels

//UPPERLOWER_FLAG - flag indicating if the hit happened in the upper (1) or lower (-1) levels







std::tuple<double,bool> hitting_times_bisection(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA)

{

double year2hour = 12*21*24;//Rescale R to hours

R = R/year2hour; // rescale R to hours

SIGMA = SIGMA/sqrt(year2hour); // rescale SIGMA to hours


/*upper/lower generation using gbm
double alpha = (1/SIGMA)*log(S_UPPER/S_CURRENT); // defines the equivalent level of the upper barrier of the GBM to the standard BM

double beta = (1/SIGMA)*log(S_LOWER/S_CURRENT); // defines the equivalent level of the upper barrier of the GBM to the standard BM

double nu = (R/SIGMA) - (SIGMA/2); // defines parameter nu


//Generate two samples of the uniform distribution

double U1 = .99*rand1();//U1 is used for the upper hit

double U2 = .99*rand1();//U2 is used for the lower hit


double inf = 1.7e308;
//Call the bisection method for both upper and lower barriers

double t_upper = root1d::secant_bisect([=](double t){return U1 - normcdf((-alpha + nu*t)/sqrt(t)) - exp(2*nu*alpha)*normcdf((-nu*t - alpha)/sqrt(t));},0.,1000.,1e-11,100,0.,inf);

double t_lower = root1d::secant_bisect([=](double t){return U2 - normcdf((beta - nu*t)/sqrt(t)) - exp(2*nu*beta)*normcdf((nu*t + beta)/sqrt(t));},0.,1000.,1e-11,100,0.,inf);

*/

double U1 = rand1();//U1 is used for the upper hit
double U2 = rand1();//U2 is used for the lower hit

double upp = normal_inv(1-U1/2);
double low = normal_inv(1-U2/2);

double t_upper = pow(((S_UPPER/S_CURRENT)-1)/(upp*SIGMA),2);
double t_lower = pow(((S_LOWER/S_CURRENT)-1)/(low*SIGMA),2);


//Defines the hitting time as the minium between the stoppig times of the upper and lower levels

double HIT_TIME = min(t_upper,t_lower);

//Defines the flag variable with true or false.

bool UPPERLOWER_FLAG = (t_upper >= t_lower)?true:false;

return std::make_tuple(HIT_TIME,UPPERLOWER_FLAG);

}



static AddInX xai_hitting_times_bisection(

	// Return type, C name, Excel name.

	FunctionX(XLL_LPOPERX, _T("?XLL_hitting_times_bisection"), _T("XLL_hitting_times_bisection"))

	.Arg(XLL_DOUBLEX, _T("S_CURRENT"), _T("Current stock price. (scalar) "),100)

	.Arg(XLL_DOUBLEX, _T("S_UPPER"), _T("Upper bound of limit. (scalar)"),105)

	.Arg(XLL_DOUBLEX, _T("S_LOWER"), _T("is the level of current price in stocks. (scalar) "),95)

	.Arg(XLL_DOUBLEX, _T("R"), _T("is the normal R.V w/ Mean 0 and Std 1. (scalar) "),0.01)

	.Arg(XLL_DOUBLEX, _T("SIGMA"), _T("is the level of volatility of the time step (dt). (scalar) "),.25)

//	.Arg(...) // add more args here

    .Volatile()

	.Category(CATEGORY)

	.FunctionHelp(_T("Generates the change in stock price, according to BlackScholes PDE."))

	.Documentation(_T("Optional documentation for function: None."))

);

LPOPERX WINAPI XLL_hitting_times_bisection(double S_CURRENT, double S_UPPER, double S_LOWER, double R, double SIGMA)

{

	#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.

	static OPERX o(2,1);



	// function body goes here

	try{

	std::tuple<double,bool> db = hitting_times_bisection(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA);

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

