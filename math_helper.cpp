#include <iostream>
#include <cmath>

// usuall, set xmid_guess to x2
double root_search(double x1, double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs,
                   double xmid_guess, double xacc)
{
  const int Maxnum_search = 50;
  int j;
  double tempx1,tempx2,dx,f,fmid,xmid,root;

  auto hermit = [] (double a, double b, double va, double vb,
            double dva, double dvb, double x, double rhs)
    {
      // use v(a), v(b), dv(a)/dt, dv(b)/dt to construct a cubic polynomial
      // basic function f1 satisfies: f1(a)=1, f1(b)=0, f1'(a)=0, f1'(b)=0
      double f1 = va*(2*x+b-3*a)*(x-b)*(x-b)/(b-a)/(b-a)/(b-a);
      // basic function f2 satisfies: f2(a)=0, f2(b)=1, f2'(a)=0, f2'(b)=0
      double f2 = vb*(3*b-2*x-a)*(x-a)*(x-a)/(b-a)/(b-a)/(b-a);
      // basic function f3 satisfies: f3(a)=0, f3(b)=0, f3'(a)=1, f3'(b)=0
      double f3 = dva*(x-a)*(x-b)*(x-b)/(b-a)/(b-a);
      // basic function f4 satisfies: f4(a)=0, f4(b)=0, f4'(a)=0, f4'(b)=1
      double f4 = dvb*(x-a)*(x-a)*(x-b)/(b-a)/(b-a);

      // the polynomial of v(x) - Vot_Threshold
      return f1 + f2 + f3 + f4 - rhs;
    };


  // for firing time case, fx1<0, fmid>0
  f    = hermit(x1,x2,fx1,fx2,dfx1,dfx2,x1,rhs);
  fmid = hermit(x1,x2,fx1,fx2,dfx1,dfx2,xmid_guess,rhs);
  /**************************************
    if (fabs(x2-x1)<xacc)
    {
      return x1;
    }
  ***************************************/
  if (f*fmid > 0) {
    //std::cerr << "root_search(): Failed to find root !" << std::endl;
    return NAN;
  }

  tempx1 = x1;
  tempx2 = xmid_guess;
  for (j=0; j<Maxnum_search; j++) {
    dx = tempx2 - tempx1;
    xmid = tempx1 + dx/2;
    fmid = hermit(x1,x2,fx1,fx2,dfx1,dfx2,xmid,rhs);
    if (fmid <= 0.0)
      tempx1 = xmid;
    else
      tempx2 = xmid;
    // the interval is small enough or already find the root
    if (fabs(fmid)<xacc) {
      root = xmid;
      return root;
    }
  }
  return xmid;
}
