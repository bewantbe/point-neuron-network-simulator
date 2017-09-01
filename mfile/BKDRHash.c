/* mex -largeArrayDims BKDRHash.c */
#include <mex.h>
#include <stdio.h>
#include <stdint.h>

/* https://www.byvoid.com/blog/string-hash-compare
 * http://www.partow.net/programming/hashfunctions/#BKDRHashFunction
 *
 * BKDR Hash Function
 * Comes from Brian Kernighan and Dennis Ritchie's book "The C Programming Language"
 */
uint32_t BKDRHash(char *str, size_t len)
{
  uint32_t seed = 131; /* 31 131 1313 13131 131313 etc.. */
  uint32_t hash = 0;

  while (len--) {
    hash = hash * seed + (*str++);
  }

  return (hash & 0x7FFFFFFF);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 1) {
    mexErrMsgTxt("Usage: BKDRHash(data)\n");
  }
  if (!mxIsNumeric(prhs[0]) && !mxIsLogical(prhs[0])
     && !mxIsChar(prhs[0])) {
    mexErrMsgTxt("Usage: BKDRHash(data); data can be numerical, logical or char array.\n");
  }

  uint32_t n = 0;
  if (mxIsChar(prhs[0])) {
    size_t sz = mxGetElementSize(prhs[0]) * mxGetNumberOfElements(prhs[0]);
    n = BKDRHash((char*)mxArrayToString(prhs[0]), sz);
  } if (mxIsSparse(prhs[0])) {
    size_t nc = mxGetN(prhs[0]);
    size_t n_elem = mxGetNzmax(prhs[0]);
    n = BKDRHash((char*)mxGetJc(prhs[0]), (nc+1)*sizeof(mwIndex));
    n = BKDRHash((char*)mxGetIr(prhs[0]), n_elem*sizeof(mwIndex))+ n*131;
    n = BKDRHash((char*)mxGetPr(prhs[0]), n_elem*sizeof(mwIndex))+ n*131;
  } else {
    size_t sz = mxGetElementSize(prhs[0]) * mxGetNumberOfElements(prhs[0]);
    n = BKDRHash((char*)mxGetData(prhs[0]), sz);
  }
  char hash_st[9];  /* 32bit number in hex form */
  sprintf(hash_st, "%.8X", n);
  plhs[0] = mxCreateString(hash_st);
}
