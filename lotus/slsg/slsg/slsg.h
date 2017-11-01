// header.h - header file for project
// Uncomment the following line to use features for Excel2007 and above.
//#define EXCEL12

#pragma once

#include "xll/xll.h"
#include <vector>

//#include "stocc.h"
#include "../../include/normal.h"
#include "../../include/secant.h"

#include <numeric>
#include <random>
#include <cfloat>
#define _USE_MATH_DEFINES
#include <cmath>

#define PI 3.1415926535

#ifndef CATEGORY
#define CATEGORY _T("Lotus")
#endif

typedef xll::traits<XLOPERX>::xcstr xcstr; // pointer to const string
typedef xll::traits<XLOPERX>::xword xword; // use for OPER and FP indices

//extern std::default_random_engine dre;



using namespace std;

inline default_random_engine& getdre()
{
	static default_random_engine dre;

	return dre;
}

//default_random_engine dre;

struct CallOption{
	double Expiration;
	double Strike;
};

struct LimitOrder{
	double amount;
	double level;
};

struct State{
	double time, position, price;
	LimitOrder limitOrder;
};

struct Stock{
	double S0, sigma;
};

inline double option_value(double S, double K, double tao, double sigma, double r = 0)
{
	double d1 = (log(S/K)+(r+pow(sigma,2)/2.0)*tao)/sigma/sqrt(tao);
	double d2 = d1 - sigma * sqrt(tao);
	double Nd1 = normal_cdf(d1);
	double Nd2 = normal_cdf(d2);
	return S*Nd1 - K*exp(-r*tao)*Nd2;
}

inline State nextHit(State oldState, double nextHittingTime, double epsilon, CallOption callOption)
{
	ensure (nextHittingTime > 0);
	ensure (epsilon > 0);
	ensure (callOption.Expiration > 0);
	ensure (callOption.Strike > 0);
	
	State newState;

	newState.time = nextHittingTime + oldState.time;
	newState.position = oldState.position + oldState.limitOrder.amount;
	newState.price = oldState.limitOrder.level;
	if (newState.position == 1) {
		newState.limitOrder.amount = -1;
		newState.limitOrder.level = callOption.Strike - epsilon;
	} else {
		newState.limitOrder.amount = 1;
		newState.limitOrder.level = callOption.Strike + epsilon;
	}

	return newState;
}


//calculate hitting time for GBM, results from Shreve stochastic 2, page 297, corollary 7.2.2
// a is limit order level, current_stock is current stock price, stock1 is the original starting point
inline double hitting_cdf(double GBM_Barrier, double current_stock, double sigma, double t, double r = 0 )
{
	ensure (GBM_Barrier > 0);
	ensure (current_stock > 0);
	ensure (sigma > 0);
	ensure (r >= 0);
	ensure (t >= 0);
	if( t == 0 )
		return 0;
	double barrier = log(GBM_Barrier/current_stock)/sigma;
	double drift = r/sigma-sigma/2.0;

	if (barrier < 0 ) { 
		drift *= -1;
		barrier *= -1;
	}
	double d1 = (-barrier+drift*t)/sqrt(t);
	double d2 = (-barrier-drift*t)/sqrt(t);

	return normal_cdf(d1)+exp(2*drift*barrier)*normal_cdf(d2);
}

//#include "../../include/secant.h"
inline double hitting_inv(double GBM_Barrier, double current_stock, double sigma, double p, double r = 0)
{
	ensure (GBM_Barrier > 0);
	ensure (current_stock > 0);
	ensure (sigma > 0);
	ensure (p >= 0);
	ensure (r >= 0);
	return root1d::secant_bisect([GBM_Barrier,current_stock, sigma, p, r](double t) { return hitting_cdf(GBM_Barrier,current_stock, sigma, t, r) - p; }, 0., 1.);
}

inline double hittingTime(double a, double current_stock, double sigma, double r = 0)
{
	ensure (a > 0);
	ensure (current_stock > 0);
	ensure (sigma > 0);
	ensure (r >= 0);

	double result = 0;
	double unif_rand = uniform_real_distribution<double>(0,1)(getdre());
	double mx_val = DBL_MAX;
	double cdf = hitting_cdf(a, current_stock, sigma, mx_val, r);

	if (unif_rand>cdf)
		result = DBL_MAX;
	else 
	{
		result = hitting_inv(a, current_stock, sigma, unif_rand, r);
		ensure( unif_rand - hitting_cdf(a,current_stock,sigma,result,r) < 10e-11 );
	}

	return result;
}

inline double double_hitting_cdf(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double t, double r = 0, double numSeries = 50 )
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (t >= 0);
	ensure (r >= 0);

	double l = log(upperBarrier) - log(lowerBarrier);
	double x = log(currentStock) - log(lowerBarrier);
	double mu = r - 0.5*sigma*sigma;           //change of measure
	double lambda_k;
	double lower_series_sum = 0;
	double upper_series_sum = 0;
	double cdf;

	double tol = 0;
	double eps = 1e-12;
	double temp1;
	double temp2;
	double count = 0;

	for (int k=1; k<=numSeries && (tol == 0 || fabs(tol) > eps); k++){
		lambda_k = 0.5*(mu*mu/sigma/sigma+k*k*PI*PI*sigma*sigma/l/l);
		temp1 = exp(-lambda_k*t)/lambda_k*k*PI*sin(k*PI*x/l);
		temp2 = exp(-lambda_k*t)/lambda_k*k*PI*sin(k*PI*(l-x)/l);
		lower_series_sum += temp1;
		upper_series_sum += temp2;
		tol = max(fabs(temp1),fabs(temp2));
		count++;
	}

	cdf = exp(-mu/sigma/sigma*x)*(sinh(mu/sigma/sigma*(l-x))/sinh(mu/sigma/sigma*l)-sigma*sigma/l/l*lower_series_sum)+
		exp(mu/sigma/sigma*(l-x))*(sinh(mu/sigma/sigma*x)/sinh(mu/sigma/sigma*l)-sigma*sigma/l/l*upper_series_sum);

	return cdf;
}


inline double double_hitting_inv(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double p, double r = 0, double numSeries = 50)
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (p >= 0);
	ensure (r >= 0);
	return root1d::secant_bisect([currentStock, upperBarrier, lowerBarrier, sigma, p, r, numSeries](double t) { return double_hitting_cdf(currentStock, upperBarrier, lowerBarrier, sigma, t, r, numSeries) - p; }, 0., 1.);
}

inline double doubleHittingTime(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double r=0, double numSeries=50)
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (r >= 0);

	double unif_rand = uniform_real_distribution<double>(0,1)(getdre());

	return double_hitting_inv(currentStock, upperBarrier, lowerBarrier, sigma, unif_rand, r, numSeries);
}

// exp(-sigma^2 t/2 + sigma B_t) ~= (1 + sigma)B_t
// S_t = s(1 + sigma B_t)
// T = inf { t : max S > s(1 + eps) or min S < s(1 - eps) } 
//   = inf { t : max B > eps/sigma or min B < -eps/sigma }
// P(T < t) = P(max B > a or min B < -a), a = eps/sigma 
// = P(max B > a) + P(min B < -a) - P(max B > a and min B < -a)
// = 2P(B > a) + 2P(B > a) - epsilon, for sigma sqrt(t) small.
// = 4P(B > a)  = 4(1 - N(a/sqrt(t))
inline double double_hitting_cdf2(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double t, double r = 0, double numSeries = 50 )
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (t >= 0);
	ensure (r >= 0);

	lowerBarrier = lowerBarrier; // not used
	double epsilon = (upperBarrier - currentStock)/currentStock;

	return (std::min)(1., 4*(1 - normal_cdf(epsilon/(sigma*sqrt(t))))); 
}

inline double double_hitting_inv2(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double p, double r = 0, double numSeries = 50)
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (p >= 0);
	ensure (r >= 0);

	lowerBarrier = lowerBarrier; // not used
	double epsilon = (upperBarrier - currentStock)/currentStock;
	double sqrt_t = epsilon/(sigma*normal_inv(1 - p/4));
	
	return sqrt_t*sqrt_t;
}

inline double doubleHittingTime2(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double r=0, double numSeries=50)
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (r >= 0);

	double unif_rand = uniform_real_distribution<double>(0,1)(getdre());

	return double_hitting_inv2(currentStock, upperBarrier, lowerBarrier, sigma, unif_rand, r, numSeries);
}

inline double doubleHittingTime3(double currentStock, double upperBarrier, double lowerBarrier, double sigma, double r=0, double numSeries=50)
{
	ensure (upperBarrier > lowerBarrier);
	ensure (lowerBarrier >= 0);
	ensure (currentStock > 0);
	ensure (sigma > 0);
	ensure (r >= 0);

	double unif_rand = uniform_real_distribution<double>(0,1)(getdre());
	
	if (unif_rand > 0.6) return double_hitting_inv(currentStock, upperBarrier, lowerBarrier, sigma, unif_rand, r, numSeries);

	else return double_hitting_inv2(currentStock, upperBarrier, lowerBarrier, sigma, unif_rand, r, numSeries);
}

inline double doubleHittingLevel(double currentStock, double upperBarrier, double lowerBarrier, double doubleHittingTime, double sigma, double r = 0, double numSeries = 50)
{
	
	double l = log(upperBarrier) - log(lowerBarrier);
	double x = log(currentStock) - log(lowerBarrier);
	double mu = r - 0.5*sigma*sigma;           //change of measure

	double unif_rand = uniform_real_distribution<double>(0,1)(getdre());

	double lambda_k;
	double upper_series_sum = 0;
	double Pr_T = double_hitting_cdf(currentStock, upperBarrier, lowerBarrier, sigma, doubleHittingTime);

	for (int k=1; k<numSeries; k++){
		lambda_k = 0.5*(mu*mu/sigma/sigma+k*k*PI*PI*sigma*sigma/l/l);
		upper_series_sum += exp(-lambda_k*doubleHittingTime)/lambda_k*k*PI*sin(k*PI*(l-x)/l);
	}

	double Pr_upperT = exp(mu/sigma/sigma*(l-x))*(sinh(mu/sigma/sigma*x)/sinh(mu/sigma/sigma*l)-sigma*sigma/l/l*upper_series_sum);
	double Pr_upper = Pr_upperT / Pr_T;

	if(unif_rand < Pr_upper) 
		return upperBarrier;
	else
		return lowerBarrier;

}


/*inline void account(size_t n, size_t m, const double* Xi, const double* X, double* A)
{
	A[0] = -inner_product(Xi, Xi + m, X, 0);

	for (size_t i = 1; i < n - 2; ++i) {
		A[i] = inner_product(Xi + m*(i - 1), Xi + m*i, X + m*i, 0); 
		A[i] -= inner_product(Xi + m*i, Xi + m*(i + 1), X + m*i, 0); 
	}

	A[n - 1] = inner_product(Xi + m*(n - 3), Xi + m*(n - 2), X + m*(n - 1), 0); 
}*/


template<class II>
inline std::vector<double> account(size_t n, size_t m, II Xi, II X)
{
	std::vector<double> A;

	A.push_back( -inner_product(Xi, Xi + m, X, 0) );

	for (size_t i = 1; i < n - 1; ++i) {
		A.push_back( inner_product(Xi + m*(i - 1), Xi + m*i, X + m*i, 0) ); 
		A[i] -= inner_product(Xi + m*i, Xi + m*(i + 1), X + m*i, 0); 
	}

	A.push_back( inner_product(Xi + m*(n - 2), Xi + m*(n - 1), X + m*(n - 1), 0)); 

	return A;
}



// Use Brownian bridge to calculate the final_price

inline double Brownian_Bridge( double level0, double level2, double t1, double t2) 
{
	ensure (t1 > 0);
	ensure (t2 > 0);
	double mean = level0+t1/t2*(level2-level0);
	double variance = t1*(t2-t1)/t2;
	double unif_rand2 = uniform_real_distribution<double>(0,1)(getdre());
	return normal_inv(unif_rand2)*sqrt(variance)+mean;
}

//increment price by time to maturity
inline double final_price_GBM(double level0, double level2, double t1, double t2, double sigma, double r = 0)
{
	ensure (level0 > 0);
	ensure (level2 > 0);
	ensure (t1 > 0);
	ensure (t2 > 0);
	ensure (sigma > 0);
	ensure (r >= 0 );
	double bblevel2 = ( log( level2 / level0 ) - ( r - 0.5 * pow(sigma, 2) ) * t2 ) / sigma ;
	double bblevel1 = Brownian_Bridge(0, bblevel2, t1, t2);

	return level0 * exp( ( r-0.5*pow(sigma,2) )*t1 + sigma*bblevel1);
}

/*
//increment price by time to maturity
inline double final_price(double level, double t)
{
	static default_random_engine dre;
	double unif_rand2 = uniform_real_distribution<double>(0,1)(dre);
	return level+normal_inv(unif_rand2)*sqrt(t);
}
*/

// we set tc_cost as a proportion to the initial stock price to find the when it is better to use limit order to hedge the option as it should
// decrease the turn over

double generate_simulation(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double width, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (width > 0);
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();
	LimitOrder lo = { 1, co.Strike + width };
	State s = { 0, 0, stock1.S0, lo };

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 
	//double cash = 0;
	double tc_cost = tc_pct * stock1.S0;

	// time 0
	Xi.push_back(s.position);
	cash = cash - s.position*(s.price + tc_cost);
	X.push_back(s.price);
	
	double t = 0;
	double nextHittingTime = hittingTime(s.limitOrder.level, s.price, stock1.sigma, r);
	ensure(nextHittingTime>0);
	t += nextHittingTime;
	while( t < co.Expiration )
	{		
		s = nextHit(s, nextHittingTime, width, co);
		Xi.push_back(s.position);
		cash = cash *exp(r*nextHittingTime) - (2*s.position - 1)*s.price - tc_cost;

		X.push_back(s.price);
		nextHittingTime = hittingTime(s.limitOrder.level, s.price, stock1.sigma, r);
		t += nextHittingTime;		
	}	
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (s.limitOrder.level - s.price) + s.price;
	//double finalPrice = final_price_GBM(s.price, s.limitOrder.level, co.Expiration-t+nextHittingTime, nextHittingTime, stock1.sigma, r);
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+s.position*s.price;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+s.position*finalPrice-abs(s.position)*tc_cost;
	return cash;
	//X.push_back(exp(r*(co.Expiration - t + nextHittingTime)) - abs(s.position) * tc_cost / cash);
}


double Simulate( const CallOption &option, const Stock &stock1, double width, double r = 0, double tc_pct = 0 )
{
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (width > 0);
	
	vector<double> Xi;
	vector<double> X;

	double hedge_diff = 0;
	hedge_diff = generate_simulation(Xi, X, option, stock1, width, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size()/2, 2, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	hedge_diff -=payoff;
	if (hedge_diff<-20)
	{
		double temp = 0;
	}
	return hedge_diff;
}


double delta_hedge_simulation(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double interval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (interval > 0);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - Nd1*tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	t += interval;

	double unirand ;

	while( t < co.Expiration )
	{	
		unirand = uniform_real_distribution<double>(0,1)(getdre());
		stock_price *=exp((r-pow(stock1.sigma,2)/2.0)*interval+stock1.sigma*normal_inv(unirand)*sqrt(interval));
		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);
		Xi.push_back(Nd1);
		cash = cash*exp(r*interval) - (Nd1-pos)*stock_price - abs(Nd1-pos)*tc_cost;
		pos = Nd1;

		X.push_back(stock_price);
		t += interval;		
	}	
	unirand = uniform_real_distribution<double>(0,1)(getdre());
	double finalPrice = stock_price * exp((r-pow(stock1.sigma,2)/2)*(co.Expiration-t+interval)+stock1.sigma*normal_inv(unirand)*sqrt(co.Expiration-t+interval));
		
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice-abs(pos)*tc_cost;
	return cash;
}

/*
double dh_simulation_priceinterval(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	double nextHittingTime = doubleHittingTime(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
  	ensure(nextHittingTime>0);
	t+=nextHittingTime;
	while( t < co.Expiration )
	{		
		stock_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
		
		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);
		Xi.push_back(Nd1);
		cash = cash*exp(r*nextHittingTime) - (Nd1-pos)*stock_price - tc_cost;
		pos = Nd1;

		X.push_back(stock_price);
		nextHittingTime = doubleHittingTime(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
		t += nextHittingTime;		
	}

	stock_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (stock_price - X.back()) + X.back();
	
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+pos*finalPrice-tc_cost;
	return cash;
}


double dh_Simulate_priceinterval( CallOption option, Stock stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);

	vector<double> Xi;
	vector<double> X;
	
	double hedge_pnl = 0;
	hedge_pnl = dh_simulation_priceinterval(Xi, X, option, stock1, pInterval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	double hedge_diff =hedge_pnl - payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}


double dh_simulation_deltainterval(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double dInterval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (dInterval > 0);
	ensure (dInterval < 1);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	double upperDelta = Nd1*(1+dInterval);
	double lowerDelta = Nd1*(1-dInterval);

	ensure (upperDelta < 1);
	ensure (lowerDelta > 0);

	double upperStock = co.Strike*exp(normal_inv(upperDelta)*stock1.sigma*sqrt(co.Expiration)-(r+pow(stock1.sigma,2)/2.0)*co.Expiration);
	double lowerStock = co.Strike*exp(normal_inv(lowerDelta)*stock1.sigma*sqrt(co.Expiration)-(r+pow(stock1.sigma,2)/2.0)*co.Expiration);
	double nextHittingTime = doubleHittingTime(stock_price, upperStock, lowerStock, stock1.sigma, r);
  	ensure(nextHittingTime>0);
	t+=nextHittingTime;
	while( t < co.Expiration )
	{		
		stock_price = doubleHittingLevel(stock_price, upperStock, lowerStock, nextHittingTime, stock1.sigma, r);
		
		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);
		Xi.push_back(Nd1);
		cash = cash*exp(r*nextHittingTime) - (Nd1-pos)*stock_price - tc_cost;
		pos = Nd1;

		X.push_back(stock_price);

		upperDelta = Nd1*(1+dInterval);
		lowerDelta = Nd1*(1-dInterval);
		upperStock = co.Strike*exp(normal_inv(upperDelta)*stock1.sigma*sqrt(co.Expiration-t)-(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t));
		lowerStock = co.Strike*exp(normal_inv(lowerDelta)*stock1.sigma*sqrt(co.Expiration-t)-(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t));

		nextHittingTime = doubleHittingTime(stock_price, upperStock, lowerStock, stock1.sigma, r);
		t += nextHittingTime;		
	}

	stock_price = doubleHittingLevel(stock_price, upperStock, lowerStock, nextHittingTime, stock1.sigma, r);
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (stock_price - X.back()) + X.back();
	
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+pos*finalPrice-tc_cost;
	return cash;
}


double dh_Simulate_deltainterval(CallOption option, Stock stock1, double dInterval, double r = 0, double tc_pct = 0)
{
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (dInterval > 0);
	ensure (dInterval < 1);

	vector<double> Xi;
	vector<double> X;
	
	double hedge_pnl = 0;
	hedge_pnl = dh_simulation_deltainterval(Xi, X, option, stock1, dInterval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	double hedge_diff =hedge_pnl - payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}

*/
double dh_simulation_gamma(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	double gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration);
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	double nextHittingTime = doubleHittingTime(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
  	ensure(nextHittingTime>0);
  	double hit_price;
  	double pos_new;

	t+=nextHittingTime;
	double count = 1;
	while( t < co.Expiration )
	{	

		hit_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
		pos_new = hit_price>stock_price?(Nd1+gamma*stock_price*pInterval):(Nd1-gamma*stock_price*pInterval);

		stock_price = hit_price;
		X.push_back(stock_price);
		Xi.push_back(pos_new);

		cash = cash*exp(r*nextHittingTime) - (pos_new-pos)*stock_price - tc_cost;
		pos = pos_new;

		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);

		gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration-t);

		nextHittingTime = doubleHittingTime(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
		t += nextHittingTime;	
		count++;
	}

	stock_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (stock_price - X.back()) + X.back();
	
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+pos*finalPrice-tc_cost;
	return cash;
}


double dh_Simulate_gamma( CallOption option, Stock stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	

	vector<double> Xi;
	vector<double> X;
	
	double hedge_pnl = 0;
	hedge_pnl = dh_simulation_gamma(Xi, X, option, stock1, pInterval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	double hedge_diff =hedge_pnl - payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}



double dh_simulation_gamma2(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	double gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration);
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	double nextHittingTime = doubleHittingTime2(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
  	ensure(nextHittingTime>0);
  	double hit_price;
  	double pos_new;

	t+=nextHittingTime;
	double count = 1;
	while( t < co.Expiration )
	{	

		hit_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
		pos_new = hit_price>stock_price?(Nd1+gamma*stock_price*pInterval):(Nd1-gamma*stock_price*pInterval);

		stock_price = hit_price;
		X.push_back(stock_price);
		Xi.push_back(pos_new);

		cash = cash*exp(r*nextHittingTime) - (pos_new-pos)*stock_price - tc_cost;
		pos = pos_new;

		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);

		gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration-t);

		nextHittingTime = doubleHittingTime2(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
		t += nextHittingTime;	
		count++;
	}

	stock_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (stock_price - X.back()) + X.back();
	
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+pos*finalPrice-tc_cost;
	return cash;
}


double dh_Simulate_gamma2( CallOption option, Stock stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	

	vector<double> Xi;
	vector<double> X;
	
	double hedge_pnl = 0;
	hedge_pnl = dh_simulation_gamma2(Xi, X, option, stock1, pInterval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	double hedge_diff =hedge_pnl - payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}

double dh_simulation_gamma3(vector<double>& Xi, vector<double>& X, const CallOption&  co, const Stock& stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	ensure (co.Expiration > 0);
	ensure (co.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	ensure ( r >= 0 );
	ensure (tc_pct >= 0);
	double tc_cost = tc_pct * stock1.S0;
	
	// n x 2 arrays of money market, stock
	Xi.clear();
	X.clear();

	// starting cash
	double cash = option_value(stock1.S0, co.Strike, co.Expiration, stock1.sigma, r); 

	// time 0
	double stock_price = stock1.S0;
	double d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*co.Expiration)/stock1.sigma/sqrt(co.Expiration);
	double Nd1 = normal_cdf(d1);
	double pos = Nd1;
	double gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration);
	Xi.push_back(Nd1);
	cash = cash-Nd1*stock_price - tc_cost;
	
	X.push_back(stock_price);
	
	double t = 0;
	double nextHittingTime = doubleHittingTime3(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
  	ensure(nextHittingTime>0);
  	double hit_price;
  	double pos_new;

	t+=nextHittingTime;
	double count = 1;
	while( t < co.Expiration )
	{	

		hit_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
		pos_new = hit_price>stock_price?(Nd1+gamma*stock_price*pInterval):(Nd1-gamma*stock_price*pInterval);

		stock_price = hit_price;
		X.push_back(stock_price);
		Xi.push_back(pos_new);

		cash = cash*exp(r*nextHittingTime) - (pos_new-pos)*stock_price - tc_cost;
		pos = pos_new;

		d1 = (log(stock_price/co.Strike)+(r+pow(stock1.sigma,2)/2.0)*(co.Expiration-t))/stock1.sigma/sqrt(co.Expiration-t);
		Nd1 = normal_cdf(d1);

		gamma = 1/sqrt(2.0*PI)*exp(-d1*d1/2.0)/stock_price/stock1.sigma/sqrt(co.Expiration-t);

		nextHittingTime = doubleHittingTime3(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), stock1.sigma, r);
		t += nextHittingTime;	
		count++;
	}

	stock_price = doubleHittingLevel(stock_price, stock_price*(1+pInterval), stock_price*(1-pInterval), nextHittingTime, stock1.sigma, r);
	
	double finalPrice = (co.Expiration - t +nextHittingTime)/nextHittingTime * (stock_price - X.back()) + X.back();
	
	X.push_back(finalPrice);
	//cash = cash*exp(r*(co.Expiration - t + interval))+pos*finalPrice;
	cash = cash*exp(r*(co.Expiration - t + nextHittingTime))+pos*finalPrice-tc_cost;
	return cash;
}


double dh_Simulate_gamma3( CallOption option, Stock stock1, double pInterval, double r = 0, double tc_pct = 0)
{
	
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (pInterval > 0);
	ensure (pInterval < 1);
	

	vector<double> Xi;
	vector<double> X;
	
	double hedge_pnl = 0;
	hedge_pnl = dh_simulation_gamma3(Xi, X, option, stock1, pInterval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	double hedge_diff =hedge_pnl - payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}



double delta_Simulate( CallOption option, Stock stock1, double interval, double r = 0, double tc_pct = 0)
{
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (interval > 0);

	vector<double> Xi;
	vector<double> X;
	
	double hedge_diff = 0;
	hedge_diff = delta_hedge_simulation(Xi, X, option, stock1, interval, r, tc_pct);
	//vector<double> A;	
	
	//A = account(X.size(), 1, Xi.begin(), X.begin());

	//hedge_diff = accumulate(A.begin(),A.end(),0);
	
	double payoff = max(X.back()-option.Strike,0);
	
	hedge_diff -=payoff;
	if (hedge_diff>10)
	{
		double temp = 0;
	}

	return hedge_diff;
}




double tracking_error_square(CallOption option, Stock stock1, double epsilon, double r = 0, double tc_pct = 0)
{
	ensure (option.Expiration > 0);
	ensure (option.Strike > 0);
	ensure (stock1.S0 > 0);
	ensure (stock1.sigma > 0);
	ensure (epsilon > 0);
	ensure (epsilon < 1);
	
	vector<double> Xi;
	vector<double> X;

	double hedge_pnl = 0;
	double payoff = 0;
	double tracking_error = 0;
	double sum = 0;
	for( int i = 0; i < 1000000; i++) {
		hedge_pnl = generate_simulation(Xi, X, option, stock1, epsilon*stock1.S0, r, tc_pct);
		payoff = max(X.back()-option.Strike,0);
		tracking_error = hedge_pnl - payoff;
		sum += pow(tracking_error, 2.0);
	};

	return sum;
}