#ifndef HEADER_NEURON_POPULATION
#define HEADER_NEURON_POPULATION

#include "common_header.h"
#include "single_neuron_dynamics.h"
#include "neuron_system_utils.h"

/*
 * Class NeuronPopulation* are abstractions for whole neuronal networks.
 * They are used to contain neuronal data and provide a simple interface
 *  to evolute and interact the network or a specific neuron.
 * It is to provide a uniform interface for simulators, bridging
 *  the single neuron models.
 */

class NeuronPopulationBase
{
public:
  // Evolution without any input. For single neuron.
  virtual void NoInteractDt(int neuron_id, double dt, double t_local, TySpikeEventVec &spike_events) = 0;
  // Evolution without any input. For population.
  virtual void NoInteractDt(double dt, double t_local, TySpikeEventVec &spike_events) {}

  // Perform an interaction
  virtual void SynapticInteraction(const TySpikeEvent &se) = 0;
  virtual void SynapticInteraction(int neuron_id, const int spike_from_id) = 0;
  virtual void InjectPoissonE(int neuron_id) = 0;
  virtual void InjectDeltaInput(int neuron_id, double strength) = 0;
  virtual void ForceReset(int neuron_id) = 0;
  virtual void DisableThreshold() = 0;
  virtual void SetThreshold(double V_thres) = 0;
  virtual void SetRefractoryTime(double t_ref) = 0;

  virtual int n_neurons() const = 0;
  virtual void operator=(const TyNeuronalDymState &neu_dym) = 0;
  virtual void ScatterCopy(const struct TyNeuronalDymState &nd,
                           const std::vector<int> &ids) = 0;

  // "Get data" functions
  virtual       TyNeuronalDymState & GetDymState() = 0;
  virtual const TyNeuronalDymState & GetDymState() const = 0;
  virtual double GetDymState(int neuron_id, int id_dym) const = 0;
  virtual const TyNeuronalParams * GetNeuronalParamsPtr() const = 0;
  virtual const Ty_Neuron_Dym_Base * GetNeuronModel() const = 0;

  // For delayed synaptic network only
  virtual double SynapticDelay() const { return 0; }
  virtual void SetSynapticDelay(double d) { };
  virtual const SparseMat * SynapticDelayNet() const { return nullptr; }
  virtual void SetSynapticDelayNet(const SparseMat &_dn) { };

  // For setting neuronal dynamical constants
  virtual void SetRisingFallingTau(const TyNeuDymParam &dym_param) { }
  virtual void SetAllDymParam(const TyArrVals &v) { }
};

class NeuronPopulationBaseCommon:
  public NeuronPopulationBase,
  public TyNeuronalParams,
  public TyNeuronalDymState
{
public:
  NeuronPopulationBaseCommon(const TyNeuronalParams &_pm, int n_var)
    :TyNeuronalParams(_pm), TyNeuronalDymState(_pm, n_var)
  {
  }

  int n_neurons() const
  {
    return Get_n_neurons();
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

// Population: delta interact, no delay, no current input
template<class TyNeu, class NBase = NeuronPopulationBaseCommon>  // neuron model
class NeuronPopulationDeltaInteractTemplate:
    public NBase
{
public:
  using NBase::StatePtr;
  using NBase::time_in_refractory;
  using NBase::n_E;
  using NBase::n_I;
  using NBase::net;
  using NBase::scee;
  using NBase::scei;
  using NBase::scie;
  using NBase::scii;
  using NBase::arr_ps;
  using NBase::dym_vals;

  TyNeu neuron_model;

  NeuronPopulationDeltaInteractTemplate(const TyNeuronalParams &_pm)
    :NBase(_pm, TyNeu::n_var)
  {}

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events) override
  {
    dbg_printf("NoInteractDt(): neuron_id = %d\n", neuron_id);
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    QUIET_STEP_CALL_INC();
    neuron_model.NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }

  void SynapticInteraction(const TySpikeEvent &se) override
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

  // FIXME: Bug here, network weight is not used here
  void SynapticInteraction(int neuron_id, const int spike_from_id) override
  {
    if (spike_from_id < n_E) {
      // Excitatory neuron fired
      if (neuron_id < n_E) {
        dym_vals(neuron_id, TyNeu::id_gEInject) += scee;
      } else {
        dym_vals(neuron_id, TyNeu::id_gEInject) += scie;
      }
    } else {
      // Inhibitory neuron fired
      if (neuron_id < n_E) {
        dym_vals(neuron_id, TyNeu::id_gIInject) += scei;
      } else {
        dym_vals(neuron_id, TyNeu::id_gIInject) += scii;
      }
    }
  }

  void ForceReset(int neuron_id) override
  {
    neuron_model.VoltHandReset(StatePtr(neuron_id));
    time_in_refractory[neuron_id] = std::numeric_limits<double>::min();
  }

  void SetRefractoryTime(double t_ref) override
  {
    neuron_model.Set_Time_Refractory(t_ref);
  }

  void SetThreshold(double V_thres) override
  {
    neuron_model.V_threshold = V_thres;
  }
  
  void DisableThreshold() override
  {
    SetThreshold(std::numeric_limits<double>::infinity());
  }

  void InjectPoissonE(int neuron_id) override
  {
    // Add Poisson input
    dym_vals(neuron_id, TyNeu::id_gEInject) += arr_ps[neuron_id];
  }

  void InjectDeltaInput(int neuron_id, double strength) override
  {
    if (strength>=0) {
      dym_vals(neuron_id, TyNeu::id_gEInject) += strength;
    } else {
      dym_vals(neuron_id, TyNeu::id_gIInject) -= strength;
    }
  }

  const Ty_Neuron_Dym_Base * GetNeuronModel() const override
  {
    return &neuron_model;
  }

  /*
  // Test and correct the neuron voltage (e.g. for imported data).
  void SaneTestVolt()
  {
    for (int j = 0; j < pm.n_total(); j++) {
      if (neu_state.dym_vals(j, p_neuron_model->Get_id_V())
          > p_neuron_model->Get_V_threshold()) {
        p_neuron_model->VoltHandReset(neu_state.StatePtr(j));
        neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();
      }
    }
  }

  // Test if the calculation blow up
  void SaneTestState() override
  {
    for (int j = 0; j < pm.n_total(); j++) {
      if ( !(fabs(neu_state.dym_vals(j, p_neuron_model->Get_id_V()))<500) ) {  // 500mV
        cerr << "\nNon-finite state value! V = " << neu_state.dym_vals(j, p_neuron_model->Get_id_V())
          << "\n  Possible reason: time step too large.\n\n";
        throw "Non-finite state value!";
      }
    }
  }

  void PrintfState(const char *const rec_path, const struct TyNeuronalDymState &neu_state)
  {
    std::ofstream fout(rec_path, std::ios_base::app);
    fout << "\n";
    fout << "State in t = " << t << "\n";
    for (int j = 0; j < neu_state.Get_n_neurons(); j++) {
      fout << "neu[" << std::setw(2) << j << "] = ";
      fout << setiosflags(std::ios::fixed) << std::setprecision(3);
      for (int k = 0; k < neu_state.Get_n_dym_vars(); k++) {
        fout << std::setw(7) << neu_state.dym_vals(j, k);
      }
      fout << "\n";
    }
    fout << std::endl;
  }
  */

};

// A model allow to set all dynamical constants for each neuron
template<class TyNeu>
class NeuronPopulationDeltaInteractHeterogeneous
  : public NeuronPopulationDeltaInteractTemplate<TyNeu>
{
  typedef NeuronPopulationDeltaInteractTemplate<TyNeu> NBase;
  using NBase::StatePtr;
  using NBase::n_neurons;
  using NBase::time_in_refractory;

public:
  std::vector<TyNeu> neuron_model_array;
  
  NeuronPopulationDeltaInteractHeterogeneous(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractTemplate<TyNeu>(_pm),
     neuron_model_array(_pm.n_total())
  { }
  
  void SetRisingFallingTau(const TyNeuDymParam &dym_param) override
  {
    if (dym_param.rows() != n_neurons()) {
      cerr << "dym_param.rows != n_neurons : " << dym_param.rows() << " != " << n_neurons() << "\n";
      throw "Inconsistent number of neurons. ";
    }
    for (int k = 0; k < n_neurons(); k++) {
      neuron_model_array[k].tau_gE_s1 = dym_param(k, 0);
      neuron_model_array[k].tau_gE    = dym_param(k, 1);
      neuron_model_array[k].tau_gI_s1 = dym_param(k, 2);
      neuron_model_array[k].tau_gI    = dym_param(k, 3);
    }
  }

  void SetAllDymParam(const TyArrVals &v) override
  {
    /*
    int neu_sz = sizeof(TyNeu);
    cout << "SetAllDymParam: sz of model = " << neu_sz << " byte.\n";
    double *vd = (double*)neuron_model_array.data();
    size_t *vi = (size_t*)neuron_model_array.data();
    for (int k = 0; k < neu_sz / sizeof(double); k++) {
      printf("v[%d] = %.16g = %lx\n", k, vd[k], vi[k]);
    }
    */
    int offset = 8;  // offset in byte, to the start of struct data member
    int neu_sz = sizeof(TyNeu) - offset;
    if (v.size() != n_neurons() * neu_sz/sizeof(double)) {
      throw "Data length inconsistent: neuronal constants. ";
    }
    // Fill in the fields in the struct of neuron model
    for (int k = 0; k < n_neurons(); k++) {
      memcpy((char*)(neuron_model_array.data() + k) + offset,
          v.data() + k * neu_sz, neu_sz);
    }
  }
  
  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events) override
  {
    dbg_printf("NoInteractDt(): neuron_id = %d\n", neuron_id);
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    QUIET_STEP_CALL_INC();
    neuron_model_array[neuron_id].NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }
};

class NeuronPopulationBaseSine: public NeuronPopulationBaseCommon
{
public:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> TySineParamVec;
  TySineParamVec sin_par;

  NeuronPopulationBaseSine(const TyNeuronalParams &_pm, int n_var)
    :NeuronPopulationBaseCommon(_pm, n_var)
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

  void SetSineFrequency(double f)
  {
    double w = 2.0*M_PI*f;
    for (int j=0; j<n_neurons(); j++) {
      sin_par(j, 1) = w;
    }
  }
};

// Population: delta interact, no delay, sine current input
template<class TyNeu>  // neuron model
class NeuronPopulationDeltaInteractSine
  :public NeuronPopulationDeltaInteractTemplate<TyNeu, NeuronPopulationBaseSine>
{
public:
  typedef NeuronPopulationDeltaInteractTemplate<TyNeu, NeuronPopulationBaseSine> NBase;
  using NBase::n_neurons;
  using NBase::StatePtr;
  using NBase::neuron_model;
  using NBase::time_in_refractory;
  using NBase::sin_par;

  NeuronPopulationDeltaInteractSine(const TyNeuronalParams &_pm)
    : NBase(_pm)
  {}

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events) override
  {
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    double *d = &sin_par(neuron_id,0);
    QUIET_STEP_CALL_INC();
    neuron_model.NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt,
        t_local, d);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }
};

// Base class for custom current input
class NeuronPopulationBaseExtI
    :public NeuronPopulationBaseCommon
{
public:
  typedef std::vector<double> TyDVec;
  TyDVec cur_param;

  NeuronPopulationBaseExtI(const TyNeuronalParams &_pm, int n_var)
    :NeuronPopulationBaseCommon(_pm, n_var)
  {}

  void SetExtI(double d)
  {
    cur_param.resize(n_neurons());
    for (int i = 0; i < n_neurons(); i++) {
      cur_param[i] = d;
    }
  }
};

// Population: delta interact, no delay, custom current input
template<class TyNeu>  // neuron model
class NeuronPopulationDeltaInteractExtI
  :public NeuronPopulationDeltaInteractTemplate<TyNeu, NeuronPopulationBaseExtI>
{
public:
  typedef NeuronPopulationDeltaInteractTemplate<TyNeu, NeuronPopulationBaseExtI> NBase;
  using NBase::n_neurons;
  using NBase::StatePtr;
  using NBase::neuron_model;
  using NBase::time_in_refractory;
  using NBase::cur_param;

  NeuronPopulationDeltaInteractExtI(const TyNeuronalParams &_pm)
    :NBase(_pm)
  {
  }

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events) override
  {
    double spike_time_local = qNaN;
    double *dym_val = StatePtr(neuron_id);
    double d = cur_param[neuron_id]; // some value
    QUIET_STEP_CALL_INC();
    neuron_model.NextStepSingleNeuronQuiet(
        dym_val, time_in_refractory[neuron_id], spike_time_local, dt,
        t_local, d);
    if (!std::isnan(spike_time_local)) {
      spike_events.emplace_back(t_local + spike_time_local, neuron_id);
    }
  }
};

// Population: delta interact, constant delay, zero current input
template<class TyNeu>
class NeuronPopulationDeltaInteractConstantDelay
  :public NeuronPopulationDeltaInteractTemplate<TyNeu>
{
public:
  double synaptic_delay;
  double SynapticDelay() const override { return synaptic_delay; }
  void SetSynapticDelay(double d) override { synaptic_delay = d; }

  NeuronPopulationDeltaInteractConstantDelay(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractTemplate<TyNeu>(_pm)
  { }
  // Let's say, it is the simulator's responsibility to relay the delayed interaction
};

// Population: delta interact, constant delay, sine current input
template<class TyNeu>
class NeuronPopulationDeltaInteractConstantDelaySine
  :public NeuronPopulationDeltaInteractSine<TyNeu> 
{
public:
  double synaptic_delay;
  double SynapticDelay() const override { return synaptic_delay; }
  void SetSynapticDelay(double d) override { synaptic_delay = d; }

  NeuronPopulationDeltaInteractConstantDelaySine(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractSine<TyNeu>(_pm)
  { }
};

// Population: delta interact, constant delay, custom current input
template<class TyNeu>
class NeuronPopulationDeltaInteractConstantDelayExtI
  :public NeuronPopulationDeltaInteractExtI<TyNeu> 
{
public:
  double synaptic_delay;
  double SynapticDelay() const override { return synaptic_delay; }
  void SetSynapticDelay(double d) override { synaptic_delay = d; }

  NeuronPopulationDeltaInteractConstantDelayExtI(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractExtI<TyNeu>(_pm)
  { }
};

// Population: delta interact, net delay, zero current input
template<class TyNeu>
class NeuronPopulationDeltaInteractNetDelay
  :public NeuronPopulationDeltaInteractTemplate<TyNeu>
{
public:
  typedef NeuronPopulationDeltaInteractTemplate<TyNeu> NBase;
  using NBase::n_neurons;
  using NBase::net;

  SparseMat synaptic_delay_net;
  double SynapticDelay() const override
  {
    throw "Non-compatible simulator - population pair. You should call SynapticDelayNet() in simulator.";
    return 0;
  }
  void SetSynapticDelay(double d) override
  {
    // Set constant delay.
    // Note: the net should be set properly first.
    synaptic_delay_net = net;
    for (int j=0; j<n_neurons(); j++) {
      for (SparseMat::InnerIterator it(synaptic_delay_net, j); it; ++it) {
        // it.row() <- j
        it.valueRef() = d;
      }
    }
  }
  const SparseMat * SynapticDelayNet() const override
  {
    return &synaptic_delay_net;
  }
  void SetSynapticDelayNet(const SparseMat &_dn) override
  {
    synaptic_delay_net = _dn;
  }

  NeuronPopulationDeltaInteractNetDelay(const TyNeuronalParams &_pm)
    :NeuronPopulationDeltaInteractTemplate<TyNeu>(_pm)
  {
    /* Set synaptic_delay_net ? */
  }
};

#endif
