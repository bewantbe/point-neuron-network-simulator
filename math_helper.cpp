#include <iostream>
#include <cmath>

/* poly/solve_cubic.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */
int
gsl_poly_solve_cubic (double a, double b, double c,
                      double *x0, double *x1, double *x2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0)
    {
      *x0 = - a / 3 ;
      *x1 = - a / 3 ;
      *x2 = - a / 3 ;
      return 3 ;
    }
  else if (CR2 == CQ3)
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          *x0 = -2 * sqrtQ  - a / 3;
          *x1 = sqrtQ - a / 3;
          *x2 = sqrtQ - a / 3;
        }
      else
        {
          *x0 = - sqrtQ  - a / 3;
          *x1 = - sqrtQ - a / 3;
          *x2 = 2 * sqrtQ - a / 3;
        }
      return 3 ;
    }
  else if (R2 < Q3)
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double ratio = sgnR * sqrt (R2 / Q3);
      double theta = acos (ratio);
      double norm = -2 * sqrt (Q);
      *x0 = norm * cos (theta / 3) - a / 3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

      /* Sort *x0, *x1, *x2 into increasing order */

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;

      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;

          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }

      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * cbrt (fabs (R) + sqrt (R2 - Q3));
      double B = Q / A ;
      *x0 = A + B - a / 3;
      return 1;
    }
}

double cubic_hermit_real_root(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs)
{
  // normalize to find root
  fx1 -= rhs;
  fx2 -= rhs;

  // normalize to x=[0,1]
  dfx1 *= x2;
  dfx2 *= x2;

  double c[4], s0=NAN, s1=NAN, s2=NAN;
  // Coefficients for hermit interpolation:
  // c[3] x^3 + c[2] x^2 + c[1] x + c[0] = 0
  c[0] = fx1;
  c[1] = dfx1;
  c[2] = -2*dfx1 - dfx2 - 3*(fx1 - fx2);
  c[3] = dfx1 + dfx2 + 2*(fx1 - fx2);
  // TODO: check c[3] == 0
  gsl_poly_solve_cubic(c[2]/c[3], c[1]/c[3], c[0]/c[3], &s0, &s1, &s2);
  if (0 <= s0 && s0 <= 1) {
    return s0 * x2;
  }
  if (0 <= s1 && s1 <= 1) {
    return s1 * x2;
  }
  if (0 <= s2 && s2 <= 1) {
    return s2 * x2;
  }
  //fprintf(stderr, "No root in this interval\n");
  return NAN;
}

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

  // Find root by binary search and Newton's iteration
  // Assume the curve is ascending.
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
