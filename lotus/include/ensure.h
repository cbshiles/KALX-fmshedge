// ensure.h - throw an exception in release mode, break in debug mode
// Copyright (c) 2006-2013 KALX, LLC. All rights reserved. No warranty is made.
//
// #define ensure(x)
// or
// #define NENSURE
// before including to turn ensure checking off

#pragma once
#include <iosfwd>
#include <stdexcept>

#ifdef NENSURE
#define ensure(x)
#endif

#ifndef ensure

#define ENSURE_HASH_(x) #x
#define ENSURE_STRZ_(x) ENSURE_HASH_(x)
#define ENSURE_FILE "file: " __FILE__
#define ENSURE_FUNC "function: " __FUNCTION__
#define ENSURE_LINE "line: " ENSURE_STRZ_(__LINE__)
#define ENSURE_SPOT ENSURE_FILE "\n" ENSURE_LINE "\n" ENSURE_FUNC

#ifdef _DEBUG
	#ifdef _WIN32 // defined for 64 bit also
		#include <Windows.h>
		#define ensure(e) if (!(e)) { DebugBreak(); }
	#elif defined(__GNUC__)
		#define ensure(e) if (!(e)) { __builtin_trap(); }
	#endif
#else // release
	#define ensure(e) if (!(e)) { \
		throw std::runtime_error(ENSURE_SPOT "\nensure: \"" #e "\" failed");}
#endif

#endif // ensure
