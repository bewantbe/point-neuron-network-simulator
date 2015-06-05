// mex randMT19937, essentially the same function as rand().
// But compatible for both Octave and Matlab

// Compile in Matlab using GCC
// mex CXXFLAGS="\$CXXFLAGS -std=c++11" CXXOPTIMFLAGS="-O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" randMT19937.cpp

// Compile in Matlab using MSVC (windows)
// mex COMPFLAGS="$COMPFLAGS" OPTIMFLAGS="$OPTIMFLAGS /fopenmp" LINKFLAGS="$LINKFLAGS /fopenmp" randMT19937.cpp

// Compile in Octave
// CXXFLAGS="-O3 -fopenmp -std=c++11"  LDFLAGS=" -fopenmp" mkoctfile --mex randMT19937.cpp

#include <mex.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <vector>
#include <random>

#define HELP_MSG \
"Usage: see help of rand()\n"

void ShowUsage()
{
  mexErrMsgTxt(HELP_MSG);
}

// declare globally, so persist during Octave running
std::random_device rd;
std::mt19937 global_urng(rd());  // gcc4.8 simd_fast_mersenne_twister_engine  // https://isocpp.org/blog/2013/03/gcc-4.8-released
std::uniform_real_distribution<> unif_dis(0, 1);

void SeedRNG(const mxArray *prhs1)
{
  global_urng.seed(0);  // seems a bug in octave, needed for global_urng.seed(sseq) to work correctly
  if (prhs1 == NULL || mxGetNumberOfElements(prhs1) == 0
      || mxIsChar(prhs1) && mxGetNumberOfElements(prhs1) == 5
      && strncmp((const char *)mxGetData(prhs1), "reset", 5) == 0) {
    // rand ("state", "reset")
    // Seed the generator by random source
    std::vector<uint32_t> li = {rd(), rd()};  // initialize with two uint
    std::seed_seq sseq(li.begin(), li.end());
    global_urng.seed(sseq);
    return;
  }
  bool b_int32;  // is int32 type seeds?
  if (mxIsChar(prhs1)) {
    // Seed the generator by characters (possiblly the state of generator)
    // note: the matlab string is not guaranteed to be '\0' terminated
    const char *st = (const char *)mxGetData(prhs1);
    mwSize n = mxGetNumberOfElements(prhs1);
    std::string cp_st(st, 0, n);
    std::stringstream ss_rng_state;
    ss_rng_state << cp_st;
    ss_rng_state >> global_urng;
    return ;
  } else if (mxIsDouble(prhs1)) {
    b_int32 = false;
  } else if (mxIsInt32(prhs1)) {
    b_int32 = true;
  } else {
    mexErrMsgTxt("type of seed not supported.");
  }
  // rand ("state", V)
  const char *px = (char*)mxGetData(prhs1);
  mwSize n = mxGetNumberOfElements(prhs1);
  std::vector<int32_t> li;
  int32_t v;
  double dv;
  for (int i=0; i<n; i++) {
    if (b_int32) {
      v = *(int32_t*)(4*i+px);
      if (v<0) {
        mexErrMsgTxt("seed_seq should be series of integers in range 0 ~ 2^31-1");
      }
    } else {
      dv = *(double*)(8*i+px);
      if (dv<0 || dv>2147483647.0) {
        mexErrMsgTxt("seed_seq should be series of integers in range 0 ~ 2^31-1");
      }
      v = (int)dv;
    }
    li.push_back(v);
  }
  std::seed_seq sseq(li.begin(), li.end());
  global_urng.seed(sseq);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  size_t m = 0, n = 0;  // number of rows and columns

  if (nrhs == 0) {
    ShowUsage();
  }
  if (mxIsChar(prhs[0])) {
    // Set or get rand state
    // Check first argument: rand ("state", ...)
    const int len = 6;  // strlen("state")
    char str[len];
    mxGetString(prhs[0], str, len);
    if (strcmp(str, "state") != 0 && strcmp(str, "seed") != 0) {
      ShowUsage();
    }
    // V = rand ("state"), V = rand("seed")
    // Output generator state
    if (nrhs == 1) {
      std::stringstream ss_rng_state;
      ss_rng_state << global_urng;
      plhs[0] = mxCreateString(ss_rng_state.str().c_str());
      return;
    }
    if (nrhs > 2) {
      ShowUsage();
    }
    // Now nrhs == 2
    SeedRNG(prhs[1]);
    return;
  }
  // Get desired dimensions
  size_t ndim = 0;            // number of dimensions
  const mwSize *dims = NULL;  // number of elements in each dimension
  std::vector<mwSize> vec_dims;
  if (nrhs == 1) {
    if (!mxIsNumeric(prhs[0])) {
      ShowUsage();
    }
    if (mxGetNumberOfElements(prhs[0]) > 1) {
      // rand ([N M ...])
      ndim = mxGetNumberOfElements(prhs[0]);  // ndim >= 2 always
      // convert double to mwSize
      vec_dims.reserve(ndim);
      double *pr = mxGetPr(prhs[0]);
      for (int k = 0; k < ndim; k++) {
        if (pr[k] < 0) {
          mexErrMsgTxt("Negative dimension size is not supported!");
          return;
        }
        vec_dims[k] = (mwSize)pr[k];
      }
      dims = vec_dims.data();
    } else {
      // rand (N)
      if (mxGetScalar(prhs[0]) < 0) {
        mexErrMsgTxt("Negative dimension size is not supported!");
        return;
      }
      vec_dims.reserve(2);
      vec_dims[0] = vec_dims[1] = mxGetScalar(prhs[0]);
      dims = vec_dims.data();
      ndim = 2;
    }
  }
  if (ndim == 0) {
    // read dimensions from each argument
    vec_dims.reserve(nrhs);
    for (int k = 0; k < nrhs; k++) {
      if (!mxIsNumeric(prhs[k]) || mxGetNumberOfElements(prhs[k])!=1) {
        ShowUsage();
      }
      if (mxGetScalar(prhs[k]) < 0) {
        mexErrMsgTxt("Negative dimension size is not supported!");
        return;
      }
      vec_dims[k] = mxGetScalar(prhs[k]);
    }
    dims = vec_dims.data();
    ndim = nrhs;
  }
  // test if there is negative dimension size
  for (int k = 0; k < vec_dims.size(); k++) {
    if (vec_dims[k] < 0) {
      mexErrMsgTxt("Negative dimension size is not supported!");
      return;
    }
  }
  // Generate random real in [0,1)
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  double *R = mxGetPr(plhs[0]);
  size_t numel = mxGetNumberOfElements(plhs[0]);
  for (size_t j = 0; j < numel; j++) {
    R[j] = unif_dis(global_urng);
  }
}

/* Tests
randMT19937
randMT19937(-1)
randMT19937(0)
randMT19937(1)
randMT19937(2)
randMT19937(2,3,4)
randMT19937([2,3,4])
c=randMT19937('state');
randMT19937('state',c);  randMT19937(1)
randMT19937('state',c);  randMT19937(1)
randMT19937('state', 123);  randMT19937(1)
randMT19937('state', 123);  randMT19937(1)
*/
