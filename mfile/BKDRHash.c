/* mex BKDRHash */
#include <mex.h>
#include <stdio.h>
#include <stdint.h>

/* https://www.byvoid.com/blog/string-hash-compare
 * http://www.partow.net/programming/hashfunctions/#BKDRHashFunction
 *
 * BKDR Hash Function
 * Comes from Brian Kernighan and Dennis Ritchie's book "The C Programming Language"
 */
uint32_t BKDRHash(char *str)
{
  uint32_t seed = 131; /* 31 131 1313 13131 131313 etc.. */
  uint32_t hash = 0;

  while (*str) {
    hash = hash * seed + (*str++);
  }

  return (hash & 0x7FFFFFFF);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 1 || mxIsChar(prhs[0]) != 1 || mxGetM(prhs[0]) > 1) {
    mexErrMsgTxt("Usage: BKDRHash(string)\n");
  }
  uint32_t n = BKDRHash(mxArrayToString(prhs[0]));
  char hash_st[9];  /* 32bit number in hex form */
  sprintf(hash_st, "%.8X", n);
  plhs[0] = mxCreateString(hash_st);
}
