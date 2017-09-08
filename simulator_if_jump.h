#include "common_header.h"
#include "neuron_system_utils.h"
#include "single_neuron_dynamics.h"
#include "poisson_generator.h"
#include "simulator_base.h"

struct Ty_IF_jump :public Ty_Neuron_Dym_Base
{
  double V_threshold  = 1.0;
  double V_reset      = 0.0;
  double G_L          = 0.05;  // ms^-1
  double T_refractory = 0.0;
  static const int n_var = 4;
  static const int id_V = 0;  // The state is V[t]
  static const int id_t = 1;
  static const int id_r = 2;  // last spike time, save this variable?
  static const int id_V_grid = 3;  // for voltage output

  bool IsInRefractory(double *dym, double t) const
  {
    return t - dym[id_r] <= T_refractory;
  }

  void ForceSpike(double *dym, double t) const
  {
    dym[id_r] = t;
    dym[id_t] = t;
    dym[id_V] = V_reset;
  }

  double GetV(double *dym, double t) const
  { // Assume V_reset == 0
    return dym[id_V] * exp(-(t - dym[id_t])*G_L);
  }

  void EvolveTo(double *dym, double t) const
  {
    dym[id_V] *= exp(-(t - dym[id_t])*G_L);
    dym[id_t] = t;
    // Maybe check if in Refractory ?
  }

  // Only forward evolve is allowed.
  void EvolveToKick(double *dym, double t, double strength) const
  {
    // Apply the poisson event.
    if (dym[id_t] != t) {
      // evolve to this time
      EvolveTo(dym, t);
    }
    if (! IsInRefractory(dym, t)) {  // TODO: no EvolveTo() also?
      dym[id_V] += strength;
    }
  }

  // override functions
  double Get_V_threshold() const override
  { return V_threshold; }

  int Get_id_V() const override
  { return id_V_grid; }

  int Get_id_gE() const override
  { return id_V_grid; }

  int Get_id_gI() const override
  { return id_V_grid; }

  int Get_id_gEInject() const override
  { return id_V_grid; }

  int Get_id_gIInject() const override
  { return id_V_grid; }

  int Get_n_dym_vars() const override
  { return n_var; }

  void Set_Time_Refractory(double t_ref) override
  { T_refractory = t_ref; }

  void NextStepSingleNeuronQuiet(
    double *dym_val,
    double &t_in_refractory,
    double &spike_time_local,
    double dt_local) const override
  {
    dym_val[id_t] += dt_local;
    if (!IsInRefractory(dym_val, dym_val[id_t]))
      dym_val[id_V] *= exp(-dt_local*G_L);
  }

  void VoltHandReset(double *dym_val) const override
  { dym_val[id_V] = V_reset; }

  const double * Get_dym_default_val() const override
  {
    static const double dym[n_var] = {0};
    return dym;
  }
};

class IFJumpPopulation :public NeuronPopulationBaseCommon
{
  Ty_IF_jump TyN;
public:
  // This class is for store data only.
  IFJumpPopulation(const TyNeuronalParams &_pm)
    :NeuronPopulationBaseCommon(_pm, Ty_IF_jump::n_var)
  {}

  void NoInteractDt(int neuron_id, double dt, double t_local, TySpikeEventVec &spike_events) override
  {}

  void SynapticInteraction(const TySpikeEvent &se) override
  {}

  void SynapticInteraction(int neuron_id, const int spike_from_id) override
  {}

  void InjectPoissonE(int neuron_id) override
  {}

  void InjectDeltaInput(int neuron_id, double strength) override
  {}

  void ForceReset(int neuron_id) override
  {}

  void DisableThreshold() override
  {}

  void SetRefractoryTime(double t_ref) override
  {}

  const Ty_Neuron_Dym_Base * GetNeuronModel() const override
  { return &TyN; }
};

class IFJumpSimulator :public NeuronSimulatorBase
{
  double t, dt;
  TyNeuronalParams pm;
  TyPoissonTimeVec poisson_time_vec;
  TySpikeEventStrengthVec input_events;
  size_t id_input_events;

  Ty_IF_jump TyN;

public:

  IFJumpSimulator(const TyNeuronalParams &_pm, double _dt)
    :pm(_pm), id_input_events(0)
  {
    t = 0;
    dt = _dt;
    poisson_time_vec.Init(pm.arr_pr, pm.arr_ps, pm.arr_pri, pm.arr_psi, t);
  }

  // Internal functions.
  void SetInputEvents(const TySpikeEventStrengthVec &_input_events)
  {
    input_events = _input_events;
    id_input_events = 0;
  }

  void AppendInputEvents()
  {
    size_t id_tail = input_events.size();
    for (size_t j = 0; j < poisson_time_vec.size(); j++) {
      const TyPoissonTimeSeq &pseq = poisson_time_vec[j];
      for (size_t k = 0; k < pseq.size(); k++)
        input_events.emplace_back(pseq[k].time, j, pseq[k].strength);
    }
    std::sort(input_events.begin() + id_tail, input_events.end());
  }

  void GeneratePoissonEvents(double t_end)
  {
    poisson_time_vec.Init(pm.arr_pr, pm.arr_ps, pm.arr_pri, pm.arr_psi, 0);
    poisson_time_vec.FillEvents(t_end, 1);

    input_events.clear();
    id_input_events = 0;
    input_events.reserve(poisson_time_vec.EventCount() + 1);
    AppendInputEvents();
    input_events.emplace_back(Inf, 0, 0);    // sealed with +Inf
  }

  void RefillPoissonEvent(double t_until, int n_least_fill)
  {
    /*poisson_time_vec.ClearAndContinueFill(t_until, n_least_fill);*/
    /*input_events.clear();*/
    /*id_input_events = 0;*/
    /*AppendInputEvents();*/

    input_events.erase(input_events.begin(),
                       input_events.begin() + id_input_events);
    id_input_events = 0;
    size_t id_tail = input_events.size();
    for (size_t j = 0; j < poisson_time_vec.size(); j++) {
      TyPoissonTimeSeq &pseq = poisson_time_vec[j];
      while (pseq.Front().time < t_until) {
        input_events.emplace_back(pseq.Front().time, j, pseq.Front().strength);
        pseq.PopAndFill();
      }
    }
    std::sort(input_events.begin() + id_tail, input_events.end());
    poisson_time_vec.SaveIdxAndClean();
  }
  
  std::vector<bool> affected;
  std::vector<double> volt_inc;
  std::vector<int> cospike_list;

  void IF_jump_simulator(TyNeuronalDymState &state, double t_end,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    const double sc_vec[2][2] = {{pm.scee, -pm.scei},{pm.scie, -pm.scii}};

    if (affected.size() != (size_t)pm.n_total()) {
      affected.resize(pm.n_total());
      volt_inc.resize(pm.n_total());
    }
    
    if (input_events.size() == 0 || input_events.back().time <= t_end) {
      RefillPoissonEvent(t_end, 1);  // TODO: +100 ?
    }

    // Loop over poisson events.
    while (id_input_events < input_events.size() &&
        input_events[id_input_events].time <= t_end) {
      const TySpikeEventStrength &e = input_events[id_input_events++];
      double *dym_val = state.StatePtr(e.id);

//      if (6 < e.time && e.time < 200) {
//        fprintf(fdbg, "e.id = %d, %f, %f\n", e.id, e.time, e.strength);
//        fprintf(fdbg, "  dym1= %f\t%f\n", dym_val[0], dym_val[1]);
//      }

      TyN.EvolveToKick(dym_val, e.time, e.strength);

      cospike_list.clear();

      // Spiked ?
      if (dym_val[TyN.id_V] > TyN.V_threshold) {
        TyN.ForceSpike(dym_val, e.time);
        cospike_list.push_back(e.id);
        ras.emplace_back(e.time, e.id);
        vec_n_spike[e.id]++;
//        if (6 < e.time && e.time < 200) {
//          fprintf(fdbg, "~~~~ spike id = %d\n", e.id);
//        }
      }

//      if (6 < e.time && e.time < 200) {
//        fprintf(fdbg, "  dym2= %f\t%f\n", dym_val[0], dym_val[1]);
//      }

      // Loop over possible cascaded spikes.
      while (! cospike_list.empty()) {
        // TODO: O(n) -> O(pn) ?
        std::fill(affected.begin(), affected.end(), false);
        std::fill(volt_inc.begin(), volt_inc.end(), 0);
        // propagate the spikes
        for (const int id : cospike_list) {
          for (SparseMat::InnerIterator it(pm.net, id); it; ++it) {
            volt_inc[it.row()] += it.value() *
              sc_vec[it.row() >= pm.n_E][id >= pm.n_E]; // TODO: expand?
            affected[it.row()] = true;
          }
        }
        cospike_list.clear();
        // Apply the increments, note down spikes.
        for (size_t j = 0; j < affected.size(); j++) {
          dym_val = state.StatePtr(j);
          if (! affected[j] || TyN.IsInRefractory(dym_val, e.time)) continue;
          TyN.EvolveToKick(dym_val, e.time, volt_inc[j]);
          if (dym_val[TyN.id_V] > TyN.V_threshold) {
            TyN.ForceSpike(dym_val, e.time);
            cospike_list.push_back(j);
            ras.emplace_back(e.time, j);
            vec_n_spike[j]++;
//            if (6 < e.time && e.time < 200) {
//              fprintf(fdbg, "~~~~ spike id = %lu\n", j);
//            }
          }
        }
      }
    }
  }
  
  // override functions
  double GetT() const override
  { return t; }

  void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike) override
  {
    t += dt;
    poisson_time_vec.SaveIdxAndClean();
    
    // Step the network.
    IF_jump_simulator(p_neu_pop->GetDymState(), t, ras, vec_n_spike);

    // Calculate voltage.
    TyNeuronalDymState &state = p_neu_pop->GetDymState();
    for (size_t j = 0; j < poisson_time_vec.size(); j++) {
      double *dym_val = state.StatePtr(j);
      dym_val[TyN.id_V_grid] = TyN.GetV(dym_val, t);
    }
  }

  // Remember to call EventTimeStrength() after set the events
  TyPoissonTimeVec & Get_poisson_time_vec() override
  {
    return poisson_time_vec;
  }

  const TyPoissonTimeVec & Get_poisson_time_vec() const override
  {
    return poisson_time_vec;
  }
};

