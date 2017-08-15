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
    if (!(stin >> id >> time)) {
      cerr << "Bad input file:\"" << path << "\"\n";
      exit(-1);
    }
    if (id >= poisson_time_vec.size()) {
      cerr << "FillPoissonEventsFromFile(): number of neuron does not match!"
        << "read id = " << id << "  number of neuron = " << poisson_time_vec.size()
        << endl;
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

void SavePoissonInput(std::ofstream &fout, TyPoissonTimeVec &poisson_time_vec, double t_step_end, const TyArrVals &arr_pr, const TyArrVals &arr_ps)
{
  poisson_time_vec.SaveIdxAndClean();
  for (size_t j = 0; j < poisson_time_vec.size(); j++) {
    TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
    const double pr = arr_pr[j];
    const double ps = arr_ps[j];
    while (poisson_time_seq.Front().time < t_step_end) {
      fout << j << " " << poisson_time_seq.Front().time
        << " " << poisson_time_seq.Front().strength << "\n";
      poisson_time_seq.PopAndFill(pr, ps);  // Next event
    }
  }
  poisson_time_vec.RestoreIdx();
}
