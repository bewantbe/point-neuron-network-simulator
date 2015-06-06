#ifndef HEADER_COMMON_HEADER
#define HEADER_COMMON_HEADER

#include <stdio.h>  // printf is good
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
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

#endif