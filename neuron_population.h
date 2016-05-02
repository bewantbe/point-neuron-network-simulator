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
};

class NeuronPopulationDeltaInteract: public NeuronPopulationBase,
                                     public TyNeuronalParams,
                                     public TyNeuronalDymState
{
  const Ty_Neuron_Dym_Base * const p_neuron_model;

public:
  NeuronPopulationDeltaInteract(const Ty_Neuron_Dym_Base * _p_neuron_model, const TyNeuronalParams &_pm)
    :TyNeuronalParams(_pm), TyNeuronalDymState(_pm, _p_neuron_model->Get_n_dym_vars()), p_neuron_model(_p_neuron_model)
  {
  }

  int n_neurons() const
  {
    return n_total();
  }

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events)
  {
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    p_neuron_model->NextStepSingleNeuronQuiet(
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
          dym_vals(it.row(), p_neuron_model->Get_id_gEInject()) += it.value() * scee;
        } else {
          dym_vals(it.row(), p_neuron_model->Get_id_gEInject()) += it.value() * scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(net, se.id); it; ++it) {
        if (it.row() < n_E) {
          dym_vals(it.row(), p_neuron_model->Get_id_gIInject()) += it.value() * scei;
        } else {
          dym_vals(it.row(), p_neuron_model->Get_id_gIInject()) += it.value() * scii;
        }
      }
    }
  }

  void ForceReset(int neuron_id)
  {
    p_neuron_model->VoltHandReset(StatePtr(neuron_id));
    time_in_refractory[neuron_id] = std::numeric_limits<double>::min();
  }

  void InjectPoissonE(int neuron_id)
  {
    // Add Poisson input
    dym_vals(neuron_id, p_neuron_model->Get_id_gEInject()) += arr_ps[neuron_id];
  }

  void operator=(const TyNeuronalDymState &neu_dym)
  {
    GetDymState() = neu_dym;
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

  /*TyNeuronalDymState & GetNeuState()*/
  /*{*/
  /*return neu_state;*/
  /*}*/
  /*const TyNeuronalDymState & GetNeuState() const*/
  /*{*/
  /*return neu_state;*/
  /*}*/
  /*double GetNeuState(int j, int id_var) const*/
  /*{*/
  /*return neu_state.dym_vals(j, id_var);*/
  /*}*/
};

#endif
