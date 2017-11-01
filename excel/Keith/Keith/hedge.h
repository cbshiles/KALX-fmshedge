// hedge.h - hedging using futures
#pragma once
#include <tuple>
#include "xll/utility/srng.h"
#include "normal.h"

namespace hedge {

	// inf {t > 0 : min B_t > lo, max B_t < hi }, hi?
	inline std::tuple<double,bool> hitting_time(double lo, double hi)
	{
		static utility::srng rng; // should use global variable

		ensure (lo < 0);
		ensure (hi > 0);

		double x;

		x = lo/normal_inv<ooura>((1 - rng.real())/2);
		double Tlo = x*x;

		x = hi/normal_inv<ooura>((1 - rng.real())/2);
		double Thi = x*x;

		return std::make_tuple(min(Tlo, Thi), Thi < Tlo); // !!! not correct
	}

	// time limit order executes, hi?
	inline std::tuple<double, bool> next_hit(double spot, double lo, double hi)
	{
		return hitting_time(lo - spot, hi - spot);
	}




	struct state {	
		double spot_;
		double vol_;
		double delta_;

		//Option Characteristics
		double k_;//Strike of the call
		double r_;
		double t_;

		double upper_;
		double lower_;

		state(double spot, double vol, double upper, double lower,double k,double r,double t)
		: spot_(spot),vol_(vol),upper_(upper),lower_(lower),k_(k),r_(r),t_(t){
		
		}

	
		
		
		void nextState(double tol){

			double hi_ = (upper_/spot_-1)/vol_;
			double lo_ = (lower_/spot_-1)/vol_;


			std::tuple<double,bool> t1 = hitting_time(lo_,hi_);
			double t = std::get<0>(t1);
			t_ -= t;
			bool upper = std::get<1>(t1);

			spot_ = upper ? upper_ : lower_;
			//s_current = upperlower *s_upper_last+(1-upperlower)*s_lower_last;
			
			upper_ = spot_*(1+tol);
			lower_ = spot_*(1-tol);



			delta_  = 	delta_call(spot_,k_,r_,vol_, t_ - tol*tol);



			

		}


	};

	






} // namespace hedge