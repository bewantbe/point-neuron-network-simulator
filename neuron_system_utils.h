#ifndef HEADER_NEURON_SYSTEM_UTILS
#define HEADER_NEURON_SYSTEM_UTILS

#include "common_header.h"

/* 
 * Structures and types (i.e. container) for neuronal system.
 */

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMat;

// Parameters about neuronal system, Especially the interactions.
struct TyNeuronalParams
{
  int n_E, n_I;    // Number of neurons, Excitatory and Inhibitory type
  SparseMat net;   // Adjacency matrix, net(i,j) means i is affect by j
  double scee, scie, scei, scii;
  TyArrVals arr_pr;   // Poisson input rate for each neuron
  TyArrVals arr_ps;   // Poisson input strength for each neuron
  // TODO: maybe add Inhibitory Poisson input?
  //TyArrVals arr_psi;

  inline int n_total() const
  { return n_E + n_I; }

  // should use this this function to initialize neurons
  void SetNumberOfNeurons(int _n_E, int _n_I)
  {
    n_E = _n_E;
    n_I = _n_I;
    arr_pr.resize(n_total());
    arr_ps.resize(n_total());
    net.resize(n_total(), n_total());
  }

  TyNeuronalParams(int _n_E, int _n_I)
  :scee(0), scie(0), scei(0), scii(0)
  {
    SetNumberOfNeurons(_n_E, _n_I);
  }
};

// All dynamical variables (V and G etc) for the neuron system
template<typename TyNeuronModel>
struct TyNeuronalDymState
{
  // Current dynamical states of neurons
  // Use RowMajor, so state variables of each neuron are continuous on memory
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dym_vals;
  // Time resided in refractory period for each neuron
  TyArrVals time_in_refractory;

  void Zeros()
  {
    dym_vals.setZero();
    for (auto &i : time_in_refractory) i = 0;
  }

  // Pointer to the state of j-th neuron
  double *StatePtr(int j)
  {
    return dym_vals.data() + j * dym_vals.cols();
  }

  TyNeuronalDymState(const TyNeuronalParams &pm)
  {
    dym_vals.resize(pm.n_total(), TyNeuronModel::n_var);
    time_in_refractory.resize(pm.n_total());
    Zeros();
  }

  TyNeuronalDymState()
  { }
};

struct TySpikeEvent
{
  double time;
  int id;

  TySpikeEvent()
  {
    time = std::numeric_limits<double>::quiet_NaN();
    id = -1;
  }

  TySpikeEvent(double _time, int _id)
  : time(_time), id(_id)
  {}

  bool operator < (const TySpikeEvent &b) const
  { return time < b.time; }

  bool operator > (const TySpikeEvent &b) const
  { return time > b.time; }

  bool operator == (const TySpikeEvent &b) const
  { return time == b.time && id == b.id; }
};

// Can be used to sort spikes increasely.
typedef std::priority_queue<
    TySpikeEvent,
    std::vector<TySpikeEvent>,
    std::greater<TySpikeEvent> > TySpikeEventQueue;

typedef std::vector< TySpikeEvent > TySpikeEventVec;

template<typename TyNeuronModel>
void FillNeuStateFromFile(TyNeuronalDymState<TyNeuronModel> &neu_dym_stat, const char *path)
{
  std::ifstream fin(path);
  double v;
  size_t j = 0;
  neu_dym_stat.dym_vals.setZero();
  while (true) {
    for (int k = 0; fin >> v, k < neu_dym_stat.dym_vals.cols(); k++) {
      neu_dym_stat.dym_vals(j, k) = v;
    }
    if (fin.fail())
      break;
    neu_dym_stat.time_in_refractory[j] = 0;
    j++;
    if (j >= neu_dym_stat.time_in_refractory.size()) {
      cerr << "read " << j << " init data" << endl;
      break;
    }
  }
}

void FillNetFromPath(TyNeuronalParams &pm, const std::string &name_net);

#endif