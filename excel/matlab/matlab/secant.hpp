// secant.hpp - 1-d root finding using secant method
#pragma once
#include <cassert>
#include <exception>
#include <limits>

#ifndef ensure
#define ensure(x) assert(x)
#endif

namespace root1d {

	template<class F, class T>
	T secant(F f, T x0, T x1, T eps = 1e-11, size_t max_iter = 100,
		T lo = -(std::numeric_limits<T>::max)(), T hi = (std::numeric_limits<T>::max)())
	{
		T f0 = f(x0);

		if (fabs(f0) <= eps)
			return x0;

		T f1 = f(x1);

		while (fabs(f1) > eps) {
			ensure (f1 != f0);

			T x_ = x1 - f1*(x1 - x0)/(f1 - f0);

			if (x_ > hi)
				x_ = hi;
			if (x_ < lo)
				x_ = lo;

			T f_ = f(x_);

			x0 = x1;
			f0 = f1;

			x1 = x_;
			f1 = f_;

			ensure (--max_iter != 0);
		}

		return x1;
	}

	template<class F, class T>
	T bisect(F f, T f0, T x0, T f1, T x1, T eps = 1e-11)
	{
		T f_, x_;

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

	template<class F, class T>
	T secant_bisect(F f, T x0, T x1, T eps = 1e-11, size_t max_iter = 100,
		T lo = -(std::numeric_limits<T>::max)(), T hi = (std::numeric_limits<T>::max)())
	{
		T f0 = f(x0);

		if ((fabs)(f0) <= eps)
			return x0;

		T f1 = f(x1);

		while (fabs(f1) > eps) {

			if (f0 * f1 < 0)
				return bisect(f, f0, x0, f1, x1, eps);

			if(f1==f0)
				f1 = f1;
			ensure (f1 != f0);

			T x_ = x1 - f1*(x1 - x0)/(f1 - f0);

			if (x_ > hi)
				x_ = hi;
			if (x_ < lo)
				x_ = lo;

			T f_ = f(x_);

			x0 = x1;
			f0 = f1;

			x1 = x_;
			f1 = f_;

			ensure (--max_iter != 0);
		}

		return x1;
	}

}; // root1d