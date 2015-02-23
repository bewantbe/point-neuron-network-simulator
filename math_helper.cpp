#include <iostream>
#include <cmath>

// usuall, set xmid_guess to x2
double root_search(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs,
                   double xmid_guess, double xacc)
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
      return (fx1*(2*rx+1) + dfx1*rx)*((rx-1)*(rx-1))
            +(fx2*(3-2*rx) + dfx2*(rx-1))*(rx*rx);
//      return fx1 + ((fx2-fx1)*(3-2*rx)*rx
//                    +(dfx1 - (dfx1+dfx2)*rx)*(1-rx))*rx;
//      return ((c3*rx+c2)*rx+c1)*rx+c0;
    };

  auto d_hermit = [&] (double rx)
    {
      return (3*c3*rx+2*c2)*rx+c1;
    };

  // for firing time case, fx1<0, fmid>0
  if (!(hermit(0) <= 0 && hermit(xmid_guess) >= 0)) {
//  if (hermit(0) * hermit(xmid_guess) > 0) {
    //std::cerr << "root_search(): Failed to find root !" << std::endl;
    return NAN;
  }

  const int Maxnum_search = 4;  // 6 + 3,  4 + 4
  double tempx1 = 0;
  double tempx2 = xmid_guess;
  double fmid,xmid;
  for (int j=0; j<Maxnum_search; j++) {
    xmid = tempx1 + (tempx2 - tempx1)/2;
    fmid = hermit(xmid);
    if (fmid <= 0.0)
      tempx1 = xmid;
    else
      tempx2 = xmid;
    // the interval is small enough or already find the root
    if (fabs(fmid)<xacc) {
      break;
//      return xmid * x2;
    }
  }
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  xmid -= hermit(xmid) / d_hermit(xmid);
  if (xmid<0 || xmid>1) {
    fprintf(stderr, "failed in Newton's iteration! Use naive guess.\n");
    xmid = 0.5;
  }
  return xmid * x2;
}
