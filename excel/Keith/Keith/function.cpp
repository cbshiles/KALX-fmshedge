// function.cpp - Rename this file and replace this description.
#include "header.h"
#include "hedge.h"


using namespace xll;

//---------------------------------Excel Add-In [1]: Implements Excel Add-In of Black-Scholes Price------------------------------------------------

static AddInX xai_function01(
	// Return type, C name, Excel name.
	FunctionX(XLL_DOUBLEX, _T("?xll_bs_price"), _T("XLL.bs_price"))
	.Arg(XLL_DOUBLEX, _T("s"), _T("is the initial Stock Price. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("K"), _T("is the Strike of Call. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("r"), _T("is the Risk-Free Rate. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("sig"), _T("is historical volatility. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("t"), _T("is time to Maturity. (scalar)"))
//	.Arg(...) // add more args here
	.Category(CATEGORY)
	.FunctionHelp(_T("Calculates the Black-Scholes Price of a call given parameters."))
	.Documentation(_T("Optional documentation: None"))
);
double WINAPI // <- Must declare all Excel C functions as WINAPI
xll_bs_price(double s,double K,double r,double sig,double t)
{
#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	//static OPERX xResult;

	// function body goes here
	double Result = price_call(s,K,r,sig,t);
	
	return Result;
}

//---------------------------------Excel Add-In [2]: Implements Excel Add-Inn of Black-Scholes Delta------------------------------------------------


static AddInX xai_function02(
	// Return type, C name, Excel name.
	FunctionX(XLL_DOUBLEX, _T("?xll_bs_delta"), _T("XLL.bs_delta"))
	.Arg(XLL_DOUBLEX, _T("s"), _T("is a initial Stock Price. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("K"), _T("is the Strike of Call. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("r"), _T("is the Risk-Free Rate. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("sig"), _T("is historical volatility. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("t"), _T("is time to Maturity. (scalar) "))
//	.Arg(...) // add more args here
	.Category(CATEGORY)
	.FunctionHelp(_T("Calculates the Black-Scholes Delta of a call given parameters."))
	.Documentation(_T("Optional documentation: None."))
);
double WINAPI // <- Must declare all Excel C functions as WINAPI
xll_bs_delta(double s,double K,double r,double sig,double t)
{
#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	//static OPERX xResult;

	// function body goes here
	double Result = delta_call(s,K,r,sig,t);
	
	return Result;
}


//---------------------------------Excel Add-In [3]: Implements Excel Add-In of Hitting Time Generation------------------------------------------------


static AddInX xai_function03(
	// Return type, C name, Excel name.
	FunctionX(XLL_LPOPERX, _T("?xll_next_hit_time"), _T("XLL.next_hit_time"))
	.Arg(XLL_DOUBLEX, _T("spot"), _T("is current spot price. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("lo"), _T("is the level of current lower limit order. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("hi"), _T("is the level of current higher limit order. (scalar) "))
	.Volatile()		//make it volatile
//	.Arg(...) // add more args here
	.Category(CATEGORY)
	.FunctionHelp(_T("Generates the next hitting time, and whether it is upper or lower limit order being placed."))
	.Documentation(_T("Optional documentation for function: None."))
);
LPOPERX WINAPI // <- Must declare all Excel C functions as WINAPI
xll_next_hit_time(double spot, double lo, double hi)
{
#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	static OPERX o(2,1);

	//using namespace hedge;
	// function body goes here
	std::tuple<double,bool> db =hedge:: next_hit(spot, lo, hi);
	o[0] = std::get<0>(db);
	o[1] = std::get<1>(db);

	return &o;
}

//---------------------------------Excel Add-In [4] & [5]: Implements Excel Add-In regarding Stocks------------------------------------------------


static AddInX xai_function04(
	// Return type, C name, Excel name.
	FunctionX(XLL_DOUBLEX, _T("?xll_priceChange_Stock"), _T("XLL.priceChange_Stock"))
	.Arg(XLL_DOUBLEX, _T("Risk-Free rate"), _T("is risk-free rate. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("dt"), _T("is the time between hedging. (scalar)"))
	.Arg(XLL_DOUBLEX, _T("current_price"), _T("is the level of current price in stocks. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("norm_rv"), _T("is the normal R.V w/ Mean 0 and Std 1. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("Sigma_Step"), _T("is the level of volatility of the time step (dt). (scalar) "))
//	.Arg(...) // add more args here
	.Category(CATEGORY)
	.FunctionHelp(_T("Generates the change in stock price, according to BlackScholes PDE."))
	.Documentation(_T("Optional documentation for function: None."))
);
double WINAPI // <- Must declare all Excel C functions as WINAPI
xll_priceChange_Stock(double r, double dt, double current_price, double norm_rv, double Sigma_Step)
{
#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	
	// function body goes here
	
	double result = r*dt*current_price+norm_rv*current_price*Sigma_Step;

	return result;
}


//---------------------------------Excel Add-In [6] & [7]: Implements Excel Add-In regarding Bond------------------------------------------------


static AddInX xai_function06(
	// Return type, C name, Excel name.
	FunctionX(XLL_DOUBLEX, _T("?xll_price_Bond"), _T("XLL.price_Bond"))
	.Arg(XLL_DOUBLEX, _T("Risk-Free rate"), _T("is risk-free rate. (scalar) "))
	.Arg(XLL_DOUBLEX, _T("Time"), _T("is current Time. (scalar)"))
	//	.Arg(...) // add more args here
	.Category(CATEGORY)
	.FunctionHelp(_T("Generates the current Bond price, assuming continuous compounding."))
	.Documentation(_T("Optional documentation for function: None."))
);
double WINAPI // <- Must declare all Excel C functions as WINAPI
xll_price_Bond(double r, double t)
{
#pragma XLLEXPORT // <- Don't forget this or your function won't get exported.
	
	// function body goes here
	
	double result = exp(r*t);

	return result;
}

static AddIn xai_state_set(
	Function(XLL_HANDLEX, "?xll_state_set", "STATE.SET")
	.Arg(XLL_DOUBLE, "spot", "is spot")
	.Arg(XLL_DOUBLE, "vol", "is vol")
	.Arg(XLL_DOUBLE, "upper", "")
	.Arg(XLL_DOUBLE, "lower" , "")
	.Arg(XLL_DOUBLE, "k", "")
	.Arg(XLL_DOUBLE, "r", "")
	.Arg(XLL_DOUBLE, "t", "")
	.Uncalced()
	.Category(CATEGORY)
	.FunctionHelp("Construct initial state")
);
HANDLEX WINAPI
xll_state_set(double spot, double vol, double upper, double lower,double k,double r,double t)
{
#pragma XLLEXPORT
	HANDLEX h(0);

	try {
		handle<hedge::state> hstate = new hedge::state(spot, vol, upper, lower,k,r,t);

		h = hstate.get();
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());

		return 0;
	}

	return h;
}
static AddIn xai_state_next(
	FunctionX(XLL_FPX,_T("?xll_state_next"), _T("STATE.NEXT"))
	.Arg(XLL_HANDLEX, _T("Handle"), _T("is a handle to the state."))
	.Arg(XLL_DOUBLE, _T("tol"), _T("is the tolerance. "))
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Advance to next state.")
);
_FP* WINAPI
xll_state_next(HANDLEX h, double tol)
{
#pragma XLLEXPORT
	static xll::FP s(7,1);

	try {
		handle<hedge::state> hstate(h);
		ensure (h);

		hstate->nextState(tol);

		s[0] = hstate->spot_;
		// ...
		s[6] = hstate->t_;
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());

		return 0;
	}

	return s.get();
}