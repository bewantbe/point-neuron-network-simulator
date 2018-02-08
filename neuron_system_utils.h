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
  TyArrVals arr_pr;   // Poisson input rate for each neuron, E type.
  TyArrVals arr_ps;   // Poisson input strength for each neuron, E type.
  TyArrVals arr_pri;   // Poisson input rate for each neuron, I type.
  TyArrVals arr_psi;   // Poisson input strength for each neuron, I type.

  inline int n_total() const
  { return n_E + n_I; }

  // should use this this function to initialize neurons
  void SetNumberOfNeurons(int _n_E, int _n_I)
  {
    n_E = _n_E;
    n_I = _n_I;
    arr_pr.resize(n_total());  std::fill(arr_pr.begin(), arr_pr.end(), 0);
    arr_ps.resize(n_total());  std::fill(arr_ps.begin(), arr_ps.end(), 0);
    arr_pri.resize(n_total()); std::fill(arr_pri.begin(), arr_pri.end(), 0);
    arr_psi.resize(n_total()); std::fill(arr_psi.begin(), arr_psi.end(), 0);
    net.resize(n_total(), n_total());
  }

  TyNeuronalParams(int _n_E, int _n_I)
  :scee(0), scie(0), scei(0), scii(0)
  {
    SetNumberOfNeurons(_n_E, _n_I);
  }
};

// All dynamical variables (V and G etc) for the neuron system
struct TyNeuronalDymState
{
  // Current dynamical states of neurons
  // Use RowMajor, so state variables of each neuron are continuous on memory
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TyDymVals;
  TyDymVals dym_vals;
  // Time resided in refractory period for each neuron
  TyArrVals time_in_refractory;

  void Zeros()
  {
    dym_vals.setZero();
    for (auto &i : time_in_refractory) i = 0;
  }

  inline int Get_n_dym_vars() const
  {
    return dym_vals.cols();
  }
  
  inline int Get_n_neurons() const
  {
    return dym_vals.rows();
  }

  // Pointer to the state of j-th neuron
  inline double *StatePtr(int j)
  {
    return dym_vals.data() + j * Get_n_dym_vars();
  }

  inline const double *StatePtr(int j) const
  {
    return dym_vals.data() + j * Get_n_dym_vars();
  }

  TyNeuronalDymState(const TyNeuronalParams &pm, const int n_dym_vars)
  {
    dym_vals.resize(pm.n_total(), n_dym_vars);
    time_in_refractory.resize(pm.n_total());
    Zeros();
  }

  // this[ids] = nd[ids]
  void ScatterCopy(const struct TyNeuronalDymState &nd,
      const std::vector<int> &ids)
  {
    for (const int &id : ids) {
      memcpy(StatePtr(id), nd.StatePtr(id),
          sizeof(double)*nd.Get_n_dym_vars());
      time_in_refractory[id] = nd.time_in_refractory[id];
    }
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
    time = qNaN;
    id = -1;
  }

  TySpikeEvent(double _time, int _id)
  : time(_time), id(_id)
  {}

  bool operator < (const TySpikeEvent &b) const
  { return time < b.time || (time == b.time && id < b.id); }

  bool operator > (const TySpikeEvent &b) const
  { return time > b.time || (time == b.time && id > b.id); }

  bool operator == (const TySpikeEvent &b) const
  { return time == b.time && id == b.id; }
};

struct TySpikeEventStrength: public TySpikeEvent
{
  double strength;

  TySpikeEventStrength()
    :TySpikeEvent()
  {
    strength = qNaN;
  }

  TySpikeEventStrength(double _time, int _id, double _strength)
    :TySpikeEvent(_time, _id), strength(_strength)
  {}
};

// Can be used to sort spikes increasely.
typedef std::priority_queue<
    TySpikeEvent,
    std::vector<TySpikeEvent>,
    std::greater<TySpikeEvent> > TySpikeEventQueue;

typedef std::vector< TySpikeEvent > TySpikeEventVec;
typedef std::vector< TySpikeEventStrength > TySpikeEventStrengthVec;

int FillNeuStateFromFile(TyNeuronalDymState &neu_dym_stat, const char *path);

void FillNetFromPath(TyNeuronalParams &pm, const std::string &name_net,
                     bool is_sparse);
SparseMat ReadNetDelay(const std::string &dn_name, const SparseMat &net);
int ReadSpikeList(TySpikeEventVec &spike_list, const char *path);

#endif
