#include "poisson_generator.h"
#include <fstream>
#include <sstream>

static const double qNaN = std::numeric_limits<double>::quiet_NaN();

// The event's times for each neuron should be ascending.
void FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps)
{
  std::ifstream fin(path);
  size_t id;
  double time;
  double strength;
  const int buf_size = 1024;
  char buf_str[buf_size];
  for (auto &i : poisson_time_vec) {
    i.clear();
  }
  while (fin.getline(buf_str, buf_size)) {
    std::istringstream stin(buf_str);
    if (!(stin >> time >> id)) {
      cerr << "Bad input file:\"" << path << "\"\n";
      exit(-1);
    }
    if (id >= poisson_time_vec.size()) {
      cerr << "FillPoissonEventsFromFile(): number of neuron does not match!" << endl;
      exit(-8);
    }
    strength = arr_ps[id];
    if (!stin.eof() && !(stin >> strength)) {
      cerr << "Bad input file:\"" << path << "\"\n";
      exit(-1);
    }
    poisson_time_vec[id].emplace_back(time, strength);
  }
  for (auto &i : poisson_time_vec) {
    i.emplace_back(qNaN, qNaN);
  }
}
