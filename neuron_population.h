#ifndef HEADER_NEURON_POPULATION
#define HEADER_NEURON_POPULATION

#include "common_header.h"
#include "single_neuron_dynamics.h"
#include "neuron_system_utils.h"

class NeuronPopulationBase
{
public:
  // For single neuron
  virtual void NoInteractDt(int neuron_id, double dt, double t_local, TySpikeEventVec &spike_events) = 0;
  // For population
  /*virtual void NoInteractDt(double dt, TySpikeEventVec &ras) = 0;*/
  // Perform an interaction
  virtual void SynapticInteraction(const TySpikeEvent &se) = 0;
  virtual void InjectPoissonE(int neuron_id) = 0;
  virtual void ForceReset(int neuron_id) = 0;
  virtual void operator=(const TyNeuronalDymState &neu_dym) = 0;
  virtual int n_neurons() const = 0;
  // "Get data" functions
  virtual       TyNeuronalDymState & GetDymState() = 0;
  virtual const TyNeuronalDymState & GetDymState() const = 0;
  virtual double GetDymState(int neuron_id, int id_dym) const = 0;
  virtual const TyNeuronalParams * GetNeuronalParamsPtr() const = 0;
  virtual void ScatterCopy(const struct TyNeuronalDymState &nd,
                           const std::vector<int> &ids) = 0;
};

template<class TyNeu>  // neuron model
class NeuronPopulationDeltaInteractTemplate:
  public NeuronPopulationBase,
  public TyNeuronalParams,
  public TyNeuronalDymState
{
public:
  TyNeu neuron_model;

  NeuronPopulationDeltaInteractTemplate(const TyNeuronalParams &_pm)
    :TyNeuronalParams(_pm), TyNeuronalDymState(_pm, TyNeu::n_var)
  {
  }

  int n_neurons() const
  {
    return Get_n_neurons();
  }

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events)
  {
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    neuron_model.NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }

  MACRO_NO_INLINE void SynapticInteraction(const TySpikeEvent &se)
  {
    if (se.id < n_E) {
      // Excitatory neuron fired
      for (SparseMat::InnerIterator it(net, se.id); it; ++it) {
        if (it.row() < n_E) {
          dym_vals(it.row(), TyNeu::id_gEInject) += it.value() * scee;
        } else {
          dym_vals(it.row(), TyNeu::id_gEInject) += it.value() * scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(net, se.id); it; ++it) {
        if (it.row() < n_E) {
          dym_vals(it.row(), TyNeu::id_gIInject) += it.value() * scei;
        } else {
          dym_vals(it.row(), TyNeu::id_gIInject) += it.value() * scii;
        }
      }
    }
  }

  void ForceReset(int neuron_id)
  {
    neuron_model.VoltHandReset(StatePtr(neuron_id));
    time_in_refractory[neuron_id] = std::numeric_limits<double>::min();
  }

  void InjectPoissonE(int neuron_id)
  {
    // Add Poisson input
    dym_vals(neuron_id, TyNeu::id_gEInject) += arr_ps[neuron_id];
  }

  void operator=(const TyNeuronalDymState &neu_dym)
  {
    GetDymState() = neu_dym;
  }

  void ScatterCopy(const struct TyNeuronalDymState &nd,
                   const std::vector<int> &ids)
  {
    GetDymState().ScatterCopy(nd, ids);
  }

  TyNeuronalDymState & GetDymState()
  {
    return *(static_cast<TyNeuronalDymState*>(this));
  }
  const TyNeuronalDymState & GetDymState() const
  {
    return *(static_cast<const TyNeuronalDymState*>(this));
  }
  TyNeuronalDymState * GetDymStatePtr()
  {
    return static_cast<TyNeuronalDymState*>(this);
  }
  const TyNeuronalDymState * GetDymStatePtr() const
  {
    return static_cast<const TyNeuronalDymState*>(this);
  }
  double GetDymState(int neuron_id, int id_dym) const
  {
    return dym_vals(neuron_id, id_dym);
  }
  const TyNeuronalParams * GetNeuronalParamsPtr() const
  {
    return static_cast<const TyNeuronalParams*>(this);
  }

};

// Neuron population with sine current input and Poisson input
template<class TyNeu>  // neuron model
class NeuronPopulationDeltaInteractSine
  :public NeuronPopulationDeltaInteractTemplate<TyNeu>
{
public:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TySineParamVec;
  TySineParamVec sin_par;

  using NeuronPopulationDeltaInteractTemplate<TyNeu>::n_neurons;
  using NeuronPopulationDeltaInteractTemplate<TyNeu>::StatePtr;
  using NeuronPopulationDeltaInteractTemplate<TyNeu>::neuron_model;
  using NeuronPopulationDeltaInteractTemplate<TyNeu>::time_in_refractory;

  NeuronPopulationDeltaInteractSine(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractTemplate<TyNeu>(_pm)
  {
    sin_par.resize(n_neurons(), 4);  // 3 or 4 both ok
    for (int j=0; j<n_neurons(); j++) {
      sin_par(j, 0) = 0;
      sin_par(j, 1) = 0;
      sin_par(j, 2) = 2*M_PI / n_neurons() * j;
    }
  }

  void SetSineAmplitude(double a)
  {
    for (int j=0; j<n_neurons(); j++) {
      sin_par(j, 0) = a;
    }
  }

  void SetSineAngularFrequency(double w)
  {
    for (int j=0; j<n_neurons(); j++) {
      sin_par(j, 1) = w;
    }
  }

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events) override
  {
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    double *d = &sin_par(neuron_id,0);
    neuron_model.NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt,
        t_local, d);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }
};

#endif
