#ifndef HEADER_MATH_HELPER
#define HEADER_MATH_HELPER

double root_search(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs);

double cubic_hermit_real_root(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double rhs);

double cubic_hermit_real_peak(double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2);

#endif
