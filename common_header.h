#ifndef HEADER_COMMON_HEADER
#define HEADER_COMMON_HEADER

#ifndef DEBUG
#define NDEBUG  // disable assert() and disable checks in Eigen
#endif
#include <cassert>

#include <stdio.h>  // printf is good
#define _USE_MATH_DEFINES  // For MSVC, so M_PI is defined
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <random>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using std::cout;
using std::cerr;
using std::endl;

#ifdef NDEBUG
#define dbg_printf(...) ((void)0);
#else
int tmp_dbg_printf(const char *format, ...);
#define dbg_printf tmp_dbg_printf
#endif

static const double Inf = std::numeric_limits<double>::infinity();
static const double qNaN = std::numeric_limits<double>::quiet_NaN();

typedef std::vector<double> TyArrVals;
typedef std::vector<TyArrVals> TyMatVals;

extern std::mt19937 rand_eng;
double g_rand();

// counter for calls to NextStepSingleNeuronQuiet() or NoInteractDt().
extern size_t quiet_step_call_counter;
#define QUIET_STEP_CALL_INC() quiet_step_call_counter++
//#define QUIET_STEP_CALL_INC() ((void)0)

// For profiling
//#define MACRO_NO_INLINE __attribute__((noinline))
#define MACRO_NO_INLINE

// Macro for detecting OS
#if defined(_WIN32) || defined(__WIN32__)
#  define _WINDOWS_
#else  // Assume Linux, no Mac OS support
#  undef _WINDOWS_
#endif

#endif
