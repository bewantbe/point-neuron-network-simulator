// g++ -Wall -O2 -std=c++11 test_file_read.cpp poisson_generator.cpp

// gcc -Wall -O2 -c external_code/dtoa.c -o dtoa.o
// g++ -Wall -O2 -std=c++11 test_file_read.cpp poisson_generator.cpp dtoa.o
#include <iostream>
#include "poisson_generator.h"

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

void FillPoissonEventsFromFile2(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps);
void FillPoissonEventsFromFile3(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps);
void FillPoissonEventsFromFile4(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps);
int main()
{
  cout << "haha\n";
  TyPoissonTimeVec poisson_time_vec;
  const char * path = "test/poi.txt";
  int n = 15;
  TyArrVals arr_pr(n), arr_ps(n), arr_pri(n), arr_psi(n);
  std::fill(arr_pr.begin(), arr_pr.end(), 0);
  std::fill(arr_ps.begin(), arr_ps.end(), 0.125);
  std::fill(arr_pri.begin(), arr_pri.end(), 0);
  std::fill(arr_psi.begin(), arr_psi.end(), 0);
  poisson_time_vec.Init(arr_pr, arr_ps, arr_pri, arr_psi, 0.0);
  FillPoissonEventsFromFile3(poisson_time_vec, path, arr_ps);
  size_t sz = 0;
  for (size_t i = 0; i < poisson_time_vec.size(); i++) {
    sz += poisson_time_vec[i].size();
  }
  cout << "sz[0] = " << poisson_time_vec[0].size() << "\n";
  cout << "size = " << sz << "\n";
  return 0;
}
