// secant.hpp - 1-d roodouble finding using secandouble method
#pragma once
#include <exception>
#include <limits>
#include "../include/ensure.h"

namespace root1d {

	template<class F>
	double secant(F f, double x0, double x1, double eps = 1e-11, size_t max_iter = 100,
		double lo = -(std::numeric_limits<double>::max)(), double hi = (std::numeric_limits<double>::max)())
	{
		double f0 = f(x0);

		if (fabs(f0) <= eps)
			return x0;

		double f1 = f(x1);

		while (fabs(f1) > eps) {
			ensure (f1 != f0);

			double x_ = x1 - f1*(x1 - x0)/(f1 - f0);

			if (x_ > hi)
				x_ = hi;
			if (x_ < lo)
				x_ = lo;

			double f_ = f(x_);

			x0 = x1;
			f0 = f1;

			x1 = x_;
			f1 = f_;

			ensure (--max_iter != 0);
		}

		return x1;
	}

	template<class F>
	double bisect(F f, double f0, double x0, double f1, double x1, double eps = 1e-11)
	{
		double f_, x_;

		ensure (f1 * f0 < 0);

		f_ = f1;
		x_ = x1;

		while (fabs(f_) > eps) {
			x_ = (x0 + x1)/2;
			f_ = f(x_);
			if (f_ * f1 < 0) {
				x0 = x_;
				f0 = f_;
			}
			else {
				x1 = x_;
				f1 = f_;
			}
		}

		return x_;
	}

	template<class F>
	double secant_bisect(F f, double x0, double x1, double eps = 1e-11, size_t max_iter = 100,
		double lo = -(std::numeric_limits<double>::max)(), double hi = (std::numeric_limits<double>::max)())
	{
		double f0 = f(x0);

		if (fabs(f0) <= eps)
			return x0;

		double f1 = f(x1);

		while (fabs(f1) > eps) {

			if (f0 * f1 < 0)
				return bisect(f, f0, x0, f1, x1, eps);

			ensure (f1 != f0);

			double x_ = x1 - f1*(x1 - x0)/(f1 - f0);

			if (x_ > hi)
				x_ = hi;
			if (x_ < lo)
				x_ = lo;

			double f_ = f(x_);

			x0 = x1;
			f0 = f1;

			x1 = x_;
			f1 = f_;

			ensure (--max_iter != 0);
		}

		return x1;
	}

}; // root1d