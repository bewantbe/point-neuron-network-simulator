#include <mex.h>
#include <stdio.h>
#include <math.h>

/**
 *  lastRASEvent(ras, n)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 2) {
    mexErrMsgTxt("Usage: lastRASEvent(ras, n)\n");
  }
  double *ras = mxGetPr(prhs[0]);
  mwSize n_events = mxGetM(prhs[0]);
  mwSize n_neurons = mxGetScalar(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(n_neurons, 1, mxREAL);
  double *t_last = mxGetPr(plhs[0]);

  int j, n_cnt = 0;
  for (j = 0; j < n_neurons; j++) t_last[j] = 0./0.;  /* NaN */
  for (j = n_events-1; j >= 0; j--) {
    mwSize i = (mwSize)ras[j] - 1;
    if (i >= n_neurons || i < 0) {
      mexErrMsgTxt("number of neurons inconsistant.");
      return;
    }
    if (isnan(t_last[i])) {
      t_last[i] = ras[n_events + j];
      n_cnt++;
      if (n_cnt == n_neurons) {
        break;
      }
    }
  }
}
