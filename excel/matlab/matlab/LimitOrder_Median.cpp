// LimitOrder_Median.cpp - Rename this file and replace this description.
#include "matlab.h"

auto rand2 = std::bind( std::uniform_real_distribution<double>(0,1), std::default_random_engine());

using namespace xll;

static AddInX xai_LimitOrder_Median(
	// Return type, C name, Excel name.
	FunctionX(XLL_FPX, _T("?xll_limitOrder_Median"), _T("XLL.limitOrder_Median"))
	// Type, name, description
	.Arg(XLL_DOUBLEX, _T("S0"), _T("Initial Value of the spot price"),100)
	.Arg(XLL_DOUBLEX, _T("K"),	_T("Strike Price of the Call Option"),100)
	.Arg(XLL_DOUBLEX, _T("R"),	_T("Risk-Free Interest Rate"),0.1)
	.Arg(XLL_DOUBLEX, _T("T"),	_T("Maturity of the Call Option [in hours]"),288)
	.Arg(XLL_DOUBLEX, _T("SIGMA"), _T("Volatility of the underlying asset [black-scholes assumption]"),0.2)
	.Arg(XLL_DOUBLEX, _T("dSmin"), _T("symmetric Distance between one limit order and the current value of the spot"),0.1)
	 .Volatile()
	.Category(CATEGORY)
	.FunctionHelp(_T("Description of what function does."))
	.Documentation(_T("Optional documentation for function."))
);
_FP* WINAPI // <- Must declare all Excel C functions as WINAPI
xll_limitOrder_Median(double S0,double K, double R, double T, double SIGMA, double dSmin)
{
#pragma XLLEXPORT // Function Exportation
	static xll::FP output(7, 1);

	try{
	
	const double year2hour = 12*21*24;	//Time conversion (24hours/day,21trading days/month,12months/year:year->hour)
	double TtoExp					= T;					//Time to Expiration
	bool optionNotExpired   = 1;					//1: Option Not Expired 0: Option Expired
	double S_Current				= S0;				//Current Stock Price
	double B_Current				=1;					//Current Bond Price
	double C_Current				= blsprice(S_Current,K,R,TtoExp/year2hour,SIGMA);//Current Value of the Call Price

	double S_Upper_Current = S_Current + dSmin;		//Current "Level" of the Upper Limit Order
	double S_Lower_Current = S_Current - dSmin;		//Current "Level" of the Lower Limit Order

	double Q_S_Current = blsdelta(S_Current,K,R,TtoExp/year2hour,SIGMA);	//Current Delta (Quantity of the Stock)
	double Q_B_Current = (C_Current - S_Current*Q_S_Current)/B_Current;		//Current Quantity of the Bond

	double Q_S_Past=0;								//Previous Quantity of Stock
	double Q_B_Past=0;								//Previous Quantity of Bond

	double PNL_Current =  0;						//PNL

	int N_Upper_Hits =0;								//Number of times the Upper Limit Order
	int N_Lower_Hits =0;

	std::tuple<double,double> o1 ;				//tuple<double,double>: 1st element=  2nd element
	double Median_Hit_Time_Upper;			//Median Hit Time for the Upper Limit Order
	double Median_Hit_Time_Lower;			//Median Hit Time for the Lower Limit Order
	double T_Delta_Upper;							//
	double T_Delta_Lower;							//

	double Delta_Upper_Current ;
	double Delta_Lower_Current ;

	double Order_Upper_Current ;
	double Order_Lower_Current ;

	std::tuple<double,bool> o2;

	double hit_Time;
	bool upperlower_Flag ;

	//Step:
	while(TtoExp>0){
		
		o1 = median_hitting_time(S_Current, S_Upper_Current, S_Lower_Current, R, SIGMA);
		Median_Hit_Time_Upper = std::get<0>(o1);
		Median_Hit_Time_Lower = std::get<1>(o1);

		T_Delta_Upper = (TtoExp>Median_Hit_Time_Upper)? TtoExp-Median_Hit_Time_Upper:0.00001;
		T_Delta_Lower = (TtoExp>Median_Hit_Time_Lower)?  TtoExp-Median_Hit_Time_Lower:0.00001;

		Delta_Upper_Current = blsdelta(S_Upper_Current,K,R,T_Delta_Upper/year2hour,SIGMA);
		Delta_Lower_Current = blsdelta(S_Lower_Current,K,R,T_Delta_Lower/year2hour,SIGMA);

		Order_Upper_Current  = Delta_Upper_Current - Q_S_Current;
		Order_Lower_Current  = Delta_Lower_Current - Q_S_Current;


		o2 = hitting_times_bisection(S_Current, S_Upper_Current, S_Lower_Current, R, SIGMA);
		hit_Time			= std::get<0>(o2);
		upperlower_Flag  = std::get<1>(o2);


		N_Upper_Hits +=upperlower_Flag?1:0;
		N_Lower_Hits +=upperlower_Flag?0:1;

		TtoExp = TtoExp - hit_Time;

		optionNotExpired = (TtoExp>0);

		//Update Current Bond Price: B_Current
		B_Current = optionNotExpired? exp((R*(T-TtoExp)/year2hour)):exp(R*T/year2hour);
		//Update Current Stock Price: S_Current
		S_Current = optionNotExpired? (S_Current =upperlower_Flag? S_Upper_Current:S_Lower_Current):S_Current;
		//Update Current Call Price: C_Current	
		C_Current = blsprice(S_Current,K,R,TtoExp/year2hour*optionNotExpired,SIGMA)*optionNotExpired;

		//Retain the Previous Quantity of Stocks and Bond held 
		Q_S_Past = Q_S_Current;
		Q_B_Past = Q_B_Current;
		//Update Quantity of Stock held Now: Q_S_Current
		Q_S_Current = optionNotExpired? (Q_S_Current = upperlower_Flag? Q_S_Past+Order_Upper_Current:Q_S_Past+Order_Lower_Current):Q_S_Past;

		Q_B_Current = (Q_B_Past*B_Current + (Q_S_Past-Q_S_Current)*S_Current)/B_Current;
		PNL_Current = ((Q_S_Current*S_Current+Q_B_Current*B_Current)-C_Current)*optionNotExpired;

		S_Upper_Current = optionNotExpired? (S_Current +dSmin): S_Upper_Current;
		S_Lower_Current = optionNotExpired? (S_Current -dSmin): S_Lower_Current;

	}
		if(optionNotExpired==0){
		
			S_Current		= S_Lower_Current +rand2()*(S_Upper_Current-S_Lower_Current);
			C_Current		= max(S_Current-K,0);
			PNL_Current  = (Q_S_Current*S_Current+Q_B_Current*B_Current)-C_Current;
		}

		output[0]=S_Current;
		output[1]=C_Current;
		output[2]=Q_S_Current;
		output[3]=Q_B_Current;
		output[4]=PNL_Current;
		output[5]=N_Upper_Hits;
		output[6]=N_Lower_Hits;
}catch(const std::exception&ex){
	XLL_ERROR(ex.what());
	return 0;
}

	return output.get();
}
