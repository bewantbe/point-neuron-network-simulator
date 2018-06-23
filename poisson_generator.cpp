#include "poisson_generator.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>  // for strtod()

// The event's times for each neuron should be ascending.
size_t FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps)
{
  std::ifstream fin(path);
  if (!fin.good()) {
    cerr << "Fail to open input file: \"" << path << "\"\n";
    return 0;
  }
  size_t id;
  size_t count_events = 0;
  double time;
  double strength;
  const int buf_size = 1024;
  char buf_str[buf_size];
  poisson_time_vec.RemoveEvents();
  while (fin.getline(buf_str, buf_size)) {
    char *pstr = buf_str;
    char *endstr;
    bool bad_convert = false;
    id = strtol(pstr, &endstr, 10);
    bad_convert |= endstr == pstr;
    pstr = endstr;
    // strtod has to be slow:
    // http://www.exploringbinary.com/how-strtod-works-and-sometimes-doesnt/
    time = strtod(pstr, &endstr);
    bad_convert |= endstr == pstr;
    if (bad_convert) {
      cerr << "Bad input file: \"" << path << "\"\n";
      exit(-1);
    }
    pstr = endstr;
    strength = strtod(pstr, &endstr);
    if (id > poisson_time_vec.size() || id == 0) {
      cerr << "FillPoissonEventsFromFile(): number of neurons does not match!\n"
        << " Read id = " << id << "  number of neurons = "
        << poisson_time_vec.size() << endl;
      exit(-8);
    }
    id -= 1;  // convert to 0-based index
    if (endstr == pstr) {
      // no strength
      strength = arr_ps[id];
    }
    poisson_time_vec[id].emplace_back(time, strength);
    count_events++;
  }
  for (auto &i : poisson_time_vec) {  // seal the queue
    i.emplace_back(Inf, 0.0);
  }
  return count_events;
}

void SavePoissonInput(std::ofstream &fout, TyPoissonTimeVec &poisson_time_vec, double t_step_end)
{
  const int buf_size = 1024;
  char str_buf[buf_size];
  poisson_time_vec.SaveIdxAndClean();
  for (size_t j = 0; j < poisson_time_vec.size(); j++) {
    TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
    while (poisson_time_seq.Front().time < t_step_end) {
      int str_len = snprintf(str_buf, buf_size, "%lu\t%.17g\t%.17g\n",
        j + 1, poisson_time_seq.Front().time, poisson_time_seq.Front().strength);
      fout.write(str_buf, std::min(str_len, buf_size));
//      fout << j + 1 << "\t" << poisson_time_seq.Front().time
//        << "\n";
//      fout << j + 1 << "\t" << poisson_time_seq.Front().time
//        << "\t" << poisson_time_seq.Front().strength << "\n";
      poisson_time_seq.PopAndFill();  // Next event
    }
  }
  poisson_time_vec.RestoreIdx();
}
