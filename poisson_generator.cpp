#include "poisson_generator.h"

// The event's times for each neuron should be ascending.
void FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path)
{
  std::ifstream fin(path);
  size_t id;
  double time;
  for (auto &i : poisson_time_vec) {
    i.clear();
  }
  while (fin >> id >> time) {
    if (id >= poisson_time_vec.size()) {
      cerr << "FillPoissonEventsFromFile(): number of neuron does not match!" << endl;
      exit(-8);
    } else {
      poisson_time_vec[id].push_back(time);
    }
  }
  for (auto &i : poisson_time_vec) {
    i.push_back(NAN);
  }
}
