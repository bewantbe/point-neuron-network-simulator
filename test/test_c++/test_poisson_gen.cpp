#include <iostream>
#include <random>
#include <algorithm>
#include <cassert>
using std::cout;
using std::cerr;
using std::endl;

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

class TyPoissonSource {
  double rate, t;
  std::exponential_distribution<> *p_exp_dis;
public:
  void Init(double _rate, double t0)
  {
    rate = _rate;
    t = t0;
    p_exp_dis = new std::exponential_distribution<>(rate);
    t = t0 + (*p_exp_dis)(rand_eng);
  }

  TyPoissonSource()
  {
    Init(1, 0);
  }

  TyPoissonSource(double _rate, double t0)
  {
    Init(_rate, t0);
  }

  double NextEventTime()
  {
    return t += (*p_exp_dis)(rand_eng);
    //return t += rand_eng();
    //return t += g_rand();
    //return t += exp(g_rand());
  }

  double CurrentEventTime() const
  { return t; }
};

int main()
{
  TyPoissonSource ps(1.0, 0.0);
  size_t nit = (size_t)4e6;
  double s = 0;
  for (size_t i = 0; i < nit; i++) {
    s += ps.NextEventTime();
  }
  cout << "s = " << s << "\n";

  return 0;
}
