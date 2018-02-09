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
public:
  Ty_IF_jump TyN;

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

  void SetThreshold(double V_thres) override
  {
    TyN.V_threshold = V_thres;
  }

  void DisableThreshold() override
  {
    TyN.V_threshold = std::numeric_limits<double>::infinity();
  }

  void SetRefractoryTime(double t_ref) override
  {
    TyN.T_refractory = t_ref;
  }

  SparseMat synaptic_delay_net;
  bool is_constant_delay_net = false;

  const Ty_Neuron_Dym_Base * GetNeuronModel() const override
  { return &TyN; }

  double SynapticDelay() const override
  {
    if (is_constant_delay_net) {
      return *synaptic_delay_net.valuePtr();
    } else {
      throw "Non-compatible simulator - population pair. You should call SynapticDelayNet() in simulator.";
    }
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
    is_constant_delay_net = true;
  }

  const SparseMat * SynapticDelayNet() const override
  {
    return &synaptic_delay_net;
  }
  void SetSynapticDelayNet(const SparseMat &_dn) override
  {
    synaptic_delay_net = _dn;
    is_constant_delay_net = false;
  }
};

class IFJumpSimulator :public NeuronSimulatorBase
{
protected:
  double t, dt;
  TyNeuronalParams pm;
  TyPoissonTimeVec poisson_time_vec;
  TySpikeEventStrengthVec input_events;
  size_t id_input_events;

  Ty_IF_jump TyN;

public:

  IFJumpSimulator(const TyNeuronalParams &_pm, double _dt, double t0)
    :pm(_pm), id_input_events(0)
  {
    t = t0;
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

      TyN.EvolveToKick(dym_val, e.time, e.strength);

      cospike_list.clear();

      // Spiked ?
      if (dym_val[TyN.id_V] > TyN.V_threshold) {
        TyN.ForceSpike(dym_val, e.time);
        cospike_list.push_back(e.id);
        ras.emplace_back(e.time, e.id);
        vec_n_spike[e.id]++;
      }

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
  
  bool ty_neuron_init = false;

  void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike) override
  {
    if (~ty_neuron_init) {
      TyN = dynamic_cast<IFJumpPopulation*>(p_neu_pop) -> TyN;
      ty_neuron_init = true;
    }
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

class IFJumpDelayNetSimulator :public IFJumpSimulator
{
protected:
  /*
  struct SmallHeadExtractorOld
  {
    const TySpikeEventStrength SEInf(Inf, 0, 0);
    const TySpikeEventStrength * p_min;
    double val;

    const TySpikeEventStrengthVec &va, &vb;
    size_t ia, ib;

    SmallHeadExtractorOld(
        const TySpikeEventStrengthVec &_va, size_t _ia,
        const TySpikeEventStrengthVec &_vb, size_t _ib)
      :va(_va), vb(_vb), ia(_ia), ib(_ib)
    {
      Next();
    }

    double Next()
    {
      if (ia >= va.size()) {
        if (ib >= vb.size()) {
          p_min = SEInf;
          val = Inf;
          return val;
        } else {
          p_min = &vb[ib++];
        }
      } else {
        if (ib >= vb.size()) {
          p_min = &va[ia++];
        } else {
          if (va[ia] < vb[ib]) {
            p_min = &va[ia++];
          } else {
            p_min = &vb[ib++];
          }
        }
      }
      val = p_min->time;
      return val;
    }

    double StepBack()  // call this when data source updated
    {
    }
  };
  */
  struct SmallHeadExtractor
  {
    const TySpikeEventStrength SEInf;
    const TySpikeEventStrength *p_min;
    double val;
    int sw;        // point to last winner
    size_t iv[2];  // point to last compared pair
    const TySpikeEventStrength *pe[2];
    const TySpikeEventStrengthVec *pv[2];

    SmallHeadExtractor(
        const TySpikeEventStrengthVec &_va, size_t _ia,
        const TySpikeEventStrengthVec &_vb, size_t _ib)
    : SEInf(Inf, 0, 0), sw(0)
    {
      init(_va, _ia, _vb, _ib);
    }

    void init(
        const TySpikeEventStrengthVec &_va, size_t _ia,
        const TySpikeEventStrengthVec &_vb, size_t _ib) {
      pv[0] = &_va;
      pv[1] = &_vb;
      iv[0] = _ia;
      iv[1] = _ib;
      ReCompare();
    }

    void ReCompare()  // call this when data source updated
    {
      pe[0] = iv[0] < pv[0]->size() ? pv[0]->data() + iv[0] : &SEInf;
      pe[1] = iv[1] < pv[1]->size() ? pv[1]->data() + iv[1] : &SEInf;
      UpdateCmp();
    }

    double Next()  // NextCurrent
    {
      iv[sw] += iv[sw] < pv[sw]->size();
      pe[sw] = iv[sw] < pv[sw]->size() ? pv[sw]->data() + iv[sw] : &SEInf;
      return UpdateCmp();
    }

    // call this if you know the modification do not alter memory
    double UpdateCmp()
    {
      sw = pe[0]->time >= pe[1]->time;
      p_min = pe[sw];
      return val = p_min->time;
    }
  };

  TySpikeEventStrengthVec intra_events;
  size_t id_intra_events = 0;

public:
  IFJumpDelayNetSimulator(const TyNeuronalParams &_pm, double _dt, double t0)
    :IFJumpSimulator(_pm, _dt, t0)
  {}

  void IF_jump_simulator(NeuronPopulationBase * p_neu_pop, double t_end,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    const auto &synaptic_delay_net = *(p_neu_pop->SynapticDelayNet());
    const auto &net = p_neu_pop->GetNeuronalParamsPtr()->net;
    TyNeuronalDymState &state = p_neu_pop->GetDymState();
    TySpikeEventStrengthVec simultaneous_events;

    const double sc_vec[2][2] = {{pm.scee, -pm.scei},{pm.scie, -pm.scii}};

    if (input_events.size() == 0 || input_events.back().time <= t_end) {
      RefillPoissonEvent(t_end, 1);
    }

    // extract the least time event.
    SmallHeadExtractor extract_min_t(input_events, id_input_events,
                                     intra_events, id_intra_events);

    while (extract_min_t.val <= t_end) {
      double t_event = extract_min_t.val;
      simultaneous_events.clear();
      simultaneous_events.push_back(*extract_min_t.p_min);
      // Get all events at the same time=t_event
      while (extract_min_t.Next() == t_event) {
        if (extract_min_t.p_min->id == simultaneous_events.back().id) {
          // Merge events if they apply to the same neuron at the same time
          simultaneous_events.back().strength += extract_min_t.p_min->strength;
        } else {
          simultaneous_events.push_back(*extract_min_t.p_min);
        }
      }
      size_t vb_end_save = intra_events.size();
      // Apply these events at a time.
      for (auto &e : simultaneous_events) {
        double *dym_val = state.StatePtr(e.id);
        TyN.EvolveToKick(dym_val, t_event, e.strength);
        if (dym_val[TyN.id_V] > TyN.V_threshold) {
          // this neuron firing
          TyN.ForceSpike(dym_val, e.time);
          ras.emplace_back(e.time, e.id);
          vec_n_spike[e.id]++;
          // Note down afected neurons.
          for (SparseMat::InnerIterator it(net, e.id),
               delay_it(synaptic_delay_net, e.id);
               it; ++it, ++delay_it) {
            intra_events.emplace_back(
                t_event + delay_it.value(),
                it.row(),
                it.value() * sc_vec[it.row() >= pm.n_E][e.id >= pm.n_E]);
          }
        }
      }
      // possible optimization: when only one neuron fired at t_event and all delay are the same, then no need to sort.
      std::sort(intra_events.begin() + vb_end_save, intra_events.end());
      // If (1) new events inserted; and (2) old non-used event exist; and
      //    (3) last old event is later than first new event.
      //    then we need a event sorting.
      if (vb_end_save < intra_events.size() &&
          vb_end_save > extract_min_t.iv[1] &&
         !(intra_events[vb_end_save-1] < intra_events[vb_end_save])) {
        std::inplace_merge(intra_events.begin()+extract_min_t.iv[1],
            intra_events.begin()+vb_end_save, intra_events.end());
      }
      // if new event inserted
      extract_min_t.ReCompare();
      // TODO: reintialize extract_min_t to deal with Inf problem.
    }
    id_input_events = extract_min_t.iv[0];
    id_intra_events = extract_min_t.iv[1];
    if (extract_min_t.iv[1] > 1000) {
      // we need to clean up the intra_event array
      intra_events.erase(intra_events.begin(), intra_events.begin() + extract_min_t.iv[1]);
      id_intra_events = 0;
    }
  }

  void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike) override
  {
    if (~ty_neuron_init) {
      TyN = dynamic_cast<IFJumpPopulation*>(p_neu_pop) -> TyN;
      ty_neuron_init = true;
    }
    t += dt;
    poisson_time_vec.SaveIdxAndClean();
    
    // Step the network.
    IF_jump_simulator(p_neu_pop, t, ras, vec_n_spike);

    // Calculate voltage.
    TyNeuronalDymState &state = p_neu_pop->GetDymState();
    for (size_t j = 0; j < poisson_time_vec.size(); j++) {
      double *dym_val = state.StatePtr(j);
      dym_val[TyN.id_V_grid] = TyN.GetV(dym_val, t);
    }
  }
};
