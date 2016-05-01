#ifndef HEADER_COMMON_HEADER
#define HEADER_COMMON_HEADER

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
#define dbg_printf printf
#endif

typedef std::vector<double> TyArrVals;

extern std::mt19937 rand_eng;
double g_rand();

//#define MACRO_NO_INLINE __attribute__((noinline))
#define MACRO_NO_INLINE

#endif