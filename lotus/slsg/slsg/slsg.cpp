// function.cpp - Rename this file and replace this description.
#include <cassert>
#include <cmath>
#include <iterator>
#include "slsg.h"

using namespace xll;

// Only run for Debug builds.
#ifdef _DEBUG

void test_hitting_time(void)
{
	// ???
//	double t;

//	t = hittingTime(1, 2);
}

void test_account(void)
{
	double Xi[] = {1, 2, 3};
	double X[] = {-1, 1, 0, 2};
	vector<double> A; // array of 3 doubles

	A=account(4, 1, Xi, X);

	// Make sure A agrees with what you expect!!!
	assert (A[0] == 1);
	assert (A[1] == -1);
	assert (A[2] == 0);
	assert (A[3] == 6);


	// Test two columns of data
	double Xi2[] = {1, 2,
		            3, 4,
				    5, 6};
	double X2[] = { 20, 10,
					30, 5,
					-10, 0,
					10, 10};
	vector<double> A2;

	A2 = account(4, 2, Xi2, X2);

	assert(A2[0] == -40);
	assert(A2[1] == -70);
	assert(A2[2] == 20);
	assert(A2[3] == 110);
	// double X2[] = 
	// double A2[6];
}

void test_optionvalue(void)
{
	double S=100;
	double K = 100;
	double T = 1;
	double sigma = 0.2;
	double value = option_value(S,K,T,sigma);
	assert(abs(value-7.965567455405804)<10e-11);
}

void test_simulate(void)
{
	double strike = 100;
	double T = 1;
	double r = 0;
	double tc_pct = 0;
	CallOption option1 = {T, strike};
	Stock stock1 = {strike,0.2};
	double width1 = 1;
	double diff = Simulate(option1,stock1,width1,r,tc_pct);
}

/*
void test_delta_simulate(void)
{
	double strike = 100;
	double T = 1;
	double r = 0;
	double tc_pct = 0;
	CallOption option1 = {T, strike};
	Stock stock1 = {strike, 0.2};
	double interval = 1.0/250;
	double diff = delta_Simulate(option1,stock1,interval,r, tc_pct);
}


void test_dh_simulate_priceinterval(void)
{
	double strike = 100;
	double T = 1;
	double r = 0;
	double tc_pct = 0;
	CallOption option1 = {T, strike};
	Stock stock1 = {strike, 0.2};
	double pInterval = 1.0/10;
	double diff = dh_Simulate_priceinterval(option1,stock1,pInterval,r, tc_pct);
}

void test_dh_simulate_deltainterval(void)
{
	double strike = 100;
	double T = 1;
	double r = 0;
	double tc_pct = 0;
	CallOption option1 = {T, strike};
	Stock stock1 = {strike, 0.2};
	double dInterval = 1.0/10;
	double diff = dh_Simulate_deltainterval(option1,stock1,dInterval,r, tc_pct);
}


int xll_test_slsg(void)
{   
	test_hitting_time();
	test_account();
	test_simulate();
	test_delta_simulate();
	test_optionvalue();
	test_dh_simulate_priceinterval();
//	test_dh_simulate_deltainterval();
//	test_dh_simulate_gamma();
	return 1;
}
static Auto<Open> xao_test_slsg(xll_test_slsg);
*/
#endif // _DEBUG


static AddIn xai_hitting_time(
	FunctionX(XLL_DOUBLE, "?xll_hitting_time", "HITTING.TIME")
	.Arg(XLL_DOUBLE, "Level", "is the hitting level. ", 101)
	.Arg(XLL_DOUBLE, "Current", "is the current price level. ", 100)
	.Arg(XLL_DOUBLE, "sigma", "is the volatility.", .20)
	.Arg(XLL_DOUBLE, "rate", "is the interest rate.", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Generate random hitting time. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_hitting_time(double a, double current, double sigma, double r)
{
#pragma XLLEXPORT
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		 result = hittingTime(a,current,sigma, r);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("HITTING.TIME: unknown error");
	}

	return result;
}

static AddIn xai_double_hitting_time(
	FunctionX(XLL_DOUBLE, "?xll_double_hitting_time", "DOUBLEHITTING.TIME")
	.Arg(XLL_DOUBLE, "currentStock", "is the current price level.", 100)
	.Arg(XLL_DOUBLE, "upperBarier", "is the upper hitting barrier.", 101)
	.Arg(XLL_DOUBLE, "lowerBarrier", "is the lower hitting barrier.", 99)
	.Arg(XLL_DOUBLE, "sigma", "is the annualized volatility.", .2)
	.Arg(XLL_DOUBLE, "rate", "is the risk-free rate.", 0)
	.Arg(XLL_DOUBLE, "numSeries", "is the number of the terms in the infinite series. ", 30)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Generate random double barrier hitting time. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_double_hitting_time(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double r, double numSeries)
{
#pragma XLLEXPORT
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = doubleHittingTime(currentStock,upperBarrier,lowerBarrier,sigma,r, numSeries);		 
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("DOUBLEHITTING.TIME: unknown error");
	}
	return result;
}
static AddIn xai_double_hitting_time2(
	FunctionX(XLL_DOUBLE, "?xll_double_hitting_time2", "DOUBLEHITTING.TIME2")
	.Arg(XLL_DOUBLE, "currentStock", "is the current price level.", 100)
	.Arg(XLL_DOUBLE, "upperBarier", "is the upper hitting barrier.", 101)
	.Arg(XLL_DOUBLE, "lowerBarrier", "is the lower hitting barrier.", 99)
	.Arg(XLL_DOUBLE, "sigma", "is the annualized volatility.", .2)
	.Arg(XLL_DOUBLE, "rate", "is the risk-free rate.", 0)
	.Arg(XLL_DOUBLE, "numSeries", "is the number of the terms in the infinite series. ", 30)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Generate random double barrier hitting time. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_double_hitting_time2(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double r, double numSeries)
{
#pragma XLLEXPORT
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = doubleHittingTime2(currentStock,upperBarrier,lowerBarrier,sigma,r, numSeries);		 
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("DOUBLEHITTING.TIME: unknown error");
	}
	return result;
}

static AddIn xai_double_hitting_level(
	FunctionX(XLL_DOUBLE, "?xll_double_hitting_level", "DOUBLEHITTING.LEVEL")
	.Arg(XLL_DOUBLE, "currentStock", "is the current price level. ")
	.Arg(XLL_DOUBLE, "upperBarier", "is the upper hitting barrier. ")
	.Arg(XLL_DOUBLE, "lowerBarrier", "is the lower hitting barrier. ")
	.Arg(XLL_DOUBLE, "doublHittingTime", "is the next double barrier hitting time. ")
	.Arg(XLL_DOUBLE, "sigma", "is the annualized volatility. ")
	.Arg(XLL_DOUBLE, "r", "is the risk-free rate. ")
	.Arg(XLL_DOUBLE, "numSeries", "is the number of the terms in the infinite series. ")
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Generate random double barrier hitting level. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_double_hitting_level(double currentStock, double upperBarrier, double lowerBarrier, double doubleHittingTime, double sigma, double r, double numSeries)
{
#pragma XLLEXPORT
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = doubleHittingLevel(currentStock, upperBarrier, lowerBarrier, doubleHittingTime, sigma, r);			
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("DOUBLEHITTING.LEVEL: unknown error");
	}
	return result;
}

static AddIn xai_simulate(
	FunctionX(XLL_DOUBLE, "?xll_simulate", "DIFF.HEDGE")
	.Arg(XLL_DOUBLE, "T", "is the maturity of option. ",1.0)
	.Arg(XLL_DOUBLE, "K", "is the strike. ",100)
	.Arg(XLL_DOUBLE, "width", "is the price window to apply limit order. ",1)
	.Arg(XLL_DOUBLE, "S0", "is the start stock price. ",100)
	.Arg(XLL_DOUBLE, "sigma", "is stock return volatility. ",0.2)
	.Arg(XLL_DOUBLE, "r", "is the interest rate. ", 0.0)
	.Arg(XLL_DOUBLE, "tc_pct", "is the percentage of tc_cost in initial stock price. ", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Calculate the hedge difference. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_simulate(double T, double K, double width, double S0, double sigma, double r, double tc_pct)
{
#pragma XLLEXPORT
	CallOption o = {T,K};
	Stock stock1 = {S0, sigma};
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = Simulate(o, stock1, width, r, tc_pct);			
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("DIFF.HEDGE: unknown error");
	}
	return result;
}



static AddIn xai_delta(
	FunctionX(XLL_DOUBLE, "?xll_delta", "DELTA.HEDGE")
	.Arg(XLL_DOUBLE, "T", "is the maturity of option. ",1.0)
	.Arg(XLL_DOUBLE, "K", "is the strike. ",100)
	.Arg(XLL_DOUBLE, "S0", "is the start stock price. ",100)
	.Arg(XLL_DOUBLE, "sigma", "is the stock return volatility. ",0.2)
	.Arg(XLL_DOUBLE, "interval", "is the hedging interval. ",1.0/250)
	.Arg(XLL_DOUBLE, "r", "is the interest rate. ", 0.0)
	.Arg(XLL_DOUBLE, "tc_pct", "is the percentage of tc_cost in initial stock price. ", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Calculate the delta hedge difference. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_delta(double T, double K, double S0, double sigma, double interval, double r, double tc_pct)
{
#pragma XLLEXPORT
	CallOption o = {T,K};
	Stock stock1 = {S0,sigma};
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = delta_Simulate(o, stock1, interval, r, tc_pct);			
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("DELTA.HEDGE: unknown error");
	}
	return result;
}


static AddIn xai_gamma(
	FunctionX(XLL_DOUBLE, "?xll_gamma", "GAMMA.HEDGE")
	.Arg(XLL_DOUBLE, "T", "is the maturity of option. ",1.0)
	.Arg(XLL_DOUBLE, "K", "is the strike. ",100)
	.Arg(XLL_DOUBLE, "S0", "is the start stock price. ",100)
	.Arg(XLL_DOUBLE, "sigma", "is the stock return volatility. ",0.2)
	.Arg(XLL_DOUBLE, "pInterval", "is the width of price change. ", 0.01)
	.Arg(XLL_DOUBLE, "r", "is the interest rate. ", 0.0)
	.Arg(XLL_DOUBLE, "tc_pct", "is the percentage of tc_cost in initial stock price. ", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Calculate the gamma-delta hedge difference. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_gamma(double T, double K, double S0, double sigma, double pInterval, double r, double tc_pct)
{
#pragma XLLEXPORT
	CallOption o = {T,K};
	Stock stock1 = {S0,sigma};
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = dh_Simulate_gamma(o, stock1, pInterval, r, tc_pct);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("GAMMA.HEDGE: unknown error");
	}
	return result;
}


static AddIn xai_gamma2(
	FunctionX(XLL_DOUBLE, "?xll_gamma2", "GAMMA.HEDGE2")
	.Arg(XLL_DOUBLE, "T", "is the maturity of option. ",1.0)
	.Arg(XLL_DOUBLE, "K", "is the strike. ",100)
	.Arg(XLL_DOUBLE, "S0", "is the start stock price. ",100)
	.Arg(XLL_DOUBLE, "sigma", "is the stock return volatility. ",0.2)
	.Arg(XLL_DOUBLE, "pInterval", "is the width of price change. ", 0.01)
	.Arg(XLL_DOUBLE, "r", "is the interest rate. ", 0.0)
	.Arg(XLL_DOUBLE, "tc_pct", "is the percentage of tc_cost in initial stock price. ", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Calculate the gamma-delta hedge difference. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_gamma2(double T, double K, double S0, double sigma, double pInterval, double r, double tc_pct)
{
#pragma XLLEXPORT
	CallOption o = {T,K};
	Stock stock1 = {S0,sigma};
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = dh_Simulate_gamma2(o, stock1, pInterval, r, tc_pct);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("GAMMA.HEDGE2: unknown error");
	}
	return result;
}

static AddIn xai_gamma3(
	FunctionX(XLL_DOUBLE, "?xll_gamma3", "GAMMA.HEDGE3")
	.Arg(XLL_DOUBLE, "T", "is the maturity of option. ",1.0)
	.Arg(XLL_DOUBLE, "K", "is the strike. ",100)
	.Arg(XLL_DOUBLE, "S0", "is the start stock price. ",100)
	.Arg(XLL_DOUBLE, "sigma", "is the stock return volatility. ",0.2)
	.Arg(XLL_DOUBLE, "pInterval", "is the width of price change. ", 0.01)
	.Arg(XLL_DOUBLE, "r", "is the interest rate. ", 0.0)
	.Arg(XLL_DOUBLE, "tc_pct", "is the percentage of tc_cost in initial stock price. ", 0.0)
	.Volatile()
	.Category(CATEGORY)
	.FunctionHelp("Calculate the gamma-delta hedge difference. ")
	.Documentation("Documentation.")
	);
double WINAPI xll_gamma3(double T, double K, double S0, double sigma, double pInterval, double r, double tc_pct)
{
#pragma XLLEXPORT
	CallOption o = {T,K};
	Stock stock1 = {S0,sigma};
	double result(std::numeric_limits<double>::quiet_NaN());

	try {
		result = dh_Simulate_gamma3(o, stock1, pInterval, r, tc_pct);
	}
	catch (const std::exception& ex) {
		XLL_ERROR(ex.what());
	}
	catch (...) {
		XLL_ERROR("GAMMA.HEDGE3: unknown error");
	}
	return result;
}