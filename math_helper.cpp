#include <iostream>
#include <cmath>

// Search the root of hermit interpolated function that
// passes zero from below to above.
// Return NAN if no root found.
// usuall, set xmid_guess to x2
double root_search(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs,
                   double xmid_guess)
{
  // normalize to x=[0,1]
  dfx1 *= x2;
  dfx2 *= x2;
  xmid_guess /= x2;

  // normalize to find root
  fx1 -= rhs;
  fx2 -= rhs;

  double c0 = fx1;
  double c1 = dfx1;
  double c2 = -2*dfx1 - dfx2 - 3*(fx1 - fx2);
  double c3 = dfx1 + dfx2 + 2*(fx1 - fx2);
  auto hermit = [&] (double rx)
    {
//      return (fx1*(2*rx+1) + dfx1*rx)*((rx-1)*(rx-1))
//            +(fx2*(3-2*rx) + dfx2*(rx-1))*(rx*rx);
      return ((c3*rx+c2)*rx+c1)*rx+c0;
    };

  auto d_hermit = [&] (double rx)
    {
      return (3*c3*rx+2*c2)*rx+c1;
    };

  // only find root in the case of monotonically increasing
  if (!(hermit(0) <= 0 && hermit(xmid_guess) >= 0)) {
    return NAN;
  }

  // Find root by binary search and newton's iteration
  // Accuracy: 2^-6 ^2 ^2 ^2 = 3.5527e-15
  const int n_iter_binary = 6;
  double tempx1 = 0;
  double tempx2 = xmid_guess;
  double fmid,xmid;
  for (int j = 0; j < n_iter_binary; j++) {
    xmid = (tempx2 + tempx1)*0.5;   // no need to take care of overflow
    fmid = hermit(xmid);
    if (fmid <= 0.0)
      tempx1 = xmid;
    else
      tempx2 = xmid;
  }
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  if (xmid<tempx1 || tempx2<xmid) {
    fprintf(stderr, "failed in Newton's iteration! Use result of binary search.\n");
    xmid = (tempx2 + tempx1)*0.5;
  }
  return xmid * x2;
}
