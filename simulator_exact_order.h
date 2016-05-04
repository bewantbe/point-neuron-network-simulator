#ifndef HEADER_SIMULATOR_EXACT_ORDER
#define HEADER_SIMULATOR_EXACT_ORDER

#include "common_header.h"
#include "poisson_generator.h"
#include "neuron_system_utils.h"
#include "simulator_base.h"

#include "neuron_population.h"

#include <iomanip>

/**
  Solver for pulse-coupled neuron model:
    After a test step (without synaptic interaction), evolve the system
    to the time of the first spike in this delta-t step, perform the synaptic
    interaction. Then start another test step for the remaining interval,
    recursively until there is no spike in the remaining interval.
  The template parameter TyNeuronModel fully describes the single neuron dynamics.

  Time cost:
    # of calls to NextStepSingleNeuronQuiet(): (1/dt + pr + 2*p*fr)*p*T
    # of synaptic interaction: (# of edges)*fr*T.
  where the fr is mean firing rate over all neurons, p=nE+nI.
*/
class NeuronSimulatorExactSpikeOrder :public NeuronSimulatorBase
{
public:
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  NeuronSimulatorExactSpikeOrder(const TyNeuronalParams &pm, double _dt)
  {
    dt = _dt;
    t = 0;
    poisson_time_vec.Init(pm.arr_pr, t);
  }

protected:
  // Evolve all neurons without synaptic interaction
  MACRO_NO_INLINE void NextStepNoInteract(NeuronPopulationBase * p_neu_pop,
              TySpikeEventVec &spike_events, double dt)
  {
    double t_step_end = t + dt;
    dbg_printf("----- Trying t = %f .. %f -------\n", t, t_step_end);
    for (int j = 0; j < p_neu_pop->n_neurons(); j++) {
      //! tmp_neu_state.dym_vals must be Row major !
      TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
      double t_local = t;
      while (poisson_time_seq.Front() < t_step_end) {
        dbg_printf("Neuron %d receive %lu-th (%lu) Poisson input at %f\n",
                   j, poisson_time_seq.id_seq, poisson_time_seq.size(),
                   poisson_time_seq.Front());
        dbg_printf("Time from %f to %f\n", t_local, poisson_time_seq.Front());
        p_neu_pop->NoInteractDt(j, poisson_time_seq.Front() - t_local, t_local, spike_events);
        t_local = poisson_time_seq.Front();
        p_neu_pop->InjectPoissonE(j);
        poisson_time_seq.PopAndFill(p_neu_pop->GetNeuronalParamsPtr()->arr_pr[j]);
      }
      dbg_printf("Time from %f to %f\n", t_local, t_step_end);
      p_neu_pop->NoInteractDt(j, t_step_end - t_local, t_local, spike_events);
      /*dbg_printf("Neuron %d @t=%f end state %f, %f, %f\n", j, t_step_end, dym_val[0], dym_val[1], dym_val[2]);*/
    }
  }

public:
  // Evolve the whole system one dt, with interaction
  MACRO_NO_INLINE void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;
    struct TyNeuronalDymState bk_neu_state;
    struct TySpikeEvent heading_spike_event(qNaN, -1);
    TySpikeEventVec spike_events;

    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    poisson_time_vec.SaveIdxAndClean();
    while (true) {
      spike_events.clear();
      bk_neu_state = p_neu_pop->GetDymState();
      NextStepNoInteract(p_neu_pop, spike_events, t_end - t);  // Try evolve till the step end
      if (spike_events.empty()) {
        t = t_end;
        break;
      } else {
        // Find out the first spike.
        heading_spike_event = *std::min_element(spike_events.begin(), spike_events.end());
        dbg_printf("Neuron [%d] spike at t = %f\n",
                   heading_spike_event.id, heading_spike_event.time);
        // Really evolve the whole system.
        poisson_time_vec.RestoreIdx();          // replay the poisson events
        spike_events.clear();
        *p_neu_pop = bk_neu_state;
        NextStepNoInteract(p_neu_pop, spike_events, heading_spike_event.time - t);
        // Record all spike events during the time step above.
        // Ideally, there should be only one event: heading_spike_event.
        bool b_heading_spike_pushed = false;
        for (size_t i = 0; i < spike_events.size(); i ++) {
          if (spike_events[i].id == heading_spike_event.id) {
            b_heading_spike_pushed = true;
          } else {
            cerr << "Unexpected spike (spike before \"first\" spike):  ["
              << i << "] id = " << spike_events[i].id
              << " time = " << spike_events[i].time << "\n";
          }
          // Here does not do `ras.emplace_back(spike_events[i])' because
          // effectively the spike interaction are done at the same time.
          ras.emplace_back(heading_spike_event.time, spike_events[i].id);
          vec_n_spike[spike_events[i].id]++;
          p_neu_pop->SynapticInteraction(spike_events[i]);
        }
        // Force the neuron to spike, if not already.
        if (!b_heading_spike_pushed) {
          p_neu_pop->ForceReset(heading_spike_event.id);
          ras.emplace_back(heading_spike_event);
          vec_n_spike[heading_spike_event.id]++;
          p_neu_pop->SynapticInteraction(heading_spike_event);
        }
        t = heading_spike_event.time;
        poisson_time_vec.SaveIdxAndClean();
      }
    }
  }

  // Test and correct the neuron voltage (e.g. for imported data).
  void SaneTestVolt()
  {
    /*for (int j = 0; j < pm.n_total(); j++) {*/
    /*if (neu_state.dym_vals(j, p_neuron_model->Get_id_V())*/
    /*> p_neuron_model->Get_V_threshold()) {*/
    /*p_neuron_model->VoltHandReset(neu_state.StatePtr(j));*/
    /*neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();*/
    /*}*/
    /*}*/
  }

  // Test if the calculation blow up
  void SaneTestState() override
  {
    /*for (int j = 0; j < pm.n_total(); j++) {*/
    /*if ( !(fabs(neu_state.dym_vals(j, p_neuron_model->Get_id_V()))<500) ) {  // 500mV*/
    /*cerr << "\nNon-finite state value! V = " << neu_state.dym_vals(j, p_neuron_model->Get_id_V())*/
    /*<< "\n  Possible reason: time step too large.\n\n";*/
    /*throw "Non-finite state value!";*/
    /*}*/
    /*}*/
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

  TyPoissonTimeVec & Get_poisson_time_vec() override
  {
    return poisson_time_vec;
  }
  const TyPoissonTimeVec & Get_poisson_time_vec() const override
  {
    return poisson_time_vec;
  }
};

/*
  The algorithm is (almost) identical to NeuronSimulatorExactSpikeOrder.
  (if neuron_model.NextStepSingleNeuronQuiet() solves the neuron exactly,
   then they produce exactly the same results)
  Should be much faster than NeuronSimulatorExactSpikeOrder when the
  network is sparse.

  Time cost:
    # of calls to NextStepSingleNeuronQuiet(): (1/dt + pr + 2*sp*fr)*p*T
    # of synaptic interaction: (# of edges)*fr*T.
  where `fr' is mean firing rate over all neurons, p=nE+nI,
        `sp' is mean number of out edges for each neuron.
*/
class NeuronSimulatorExactSpikeOrderSparse
: public NeuronSimulatorExactSpikeOrder
{
public:

  NeuronSimulatorExactSpikeOrderSparse(const TyNeuronalParams &_pm, double _dt)
    :NeuronSimulatorExactSpikeOrder(_pm, _dt)
  {
  }

protected:
  // Evolve all neurons without synaptic interaction
  MACRO_NO_INLINE void NextStepNoInteractToTime(
      NeuronPopulationBase * p_neu_pop,
      const TyArrVals &bk_state_time,
      const std::vector<int> &ids_affected,
      TySpikeEventVec &spike_events,
      double t_step_end)
  {
    for (size_t jj = 0; jj < ids_affected.size(); jj++) {
      int j = ids_affected[jj];
      dbg_printf("----- NextStepNoInteractToTime(), Neuron %d, t = %f .. %f\n", j, bk_state_time[j], t_step_end);
      //! tmp_neu_state.dym_vals must be Row major !
      TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
      double t_local = bk_state_time[j];
      while (poisson_time_seq.Front() < t_step_end) {
        dbg_printf("  receive %lu-th (%lu) Poisson input at %f\n",
                   poisson_time_seq.id_seq, poisson_time_seq.size(),
                   poisson_time_seq.Front());
        dbg_printf("  time from %f to %f\n", t_local, poisson_time_seq.Front());
        p_neu_pop->NoInteractDt(j,
            poisson_time_seq.Front() - t_local, t_local, spike_events);
        t_local = poisson_time_seq.Front();
        p_neu_pop->InjectPoissonE(j);
        poisson_time_seq.PopAndFill(p_neu_pop->GetNeuronalParamsPtr()->arr_pr[j]);
      }
      dbg_printf("  time from %f to %f\n", t_local, t_step_end);
      p_neu_pop->NoInteractDt(j,
          t_step_end - t_local, t_local, spike_events);
      /*dbg_printf("  t=%f end state: %f, %f, %f\n", t_step_end, dym_val[0], dym_val[1], dym_val[2]);*/
    }
  }

public:

  // Evolve the whole system one dt, with interaction.
  // Save spike events in this `dt' in `ras',
  // and add spike count to `vec_n_spike'.
  // This version also roll back the firing neuron

  std::vector<int> ids_affected;             // index of affected neurons
  std::vector<bool> bool_affected;           // hash for `ids_affected'

  // `bk_neu_state' is state at time `bk_state_time'
  TyArrVals bk_state_time;
  struct TyNeuronalDymState bk_neu_state;

  MACRO_NO_INLINE void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;

    struct TySpikeEvent heading_spike_event(qNaN, -1);
    ids_affected.clear();
    bool_affected.clear();
    for (int i = 0; i < p_neu_pop->n_neurons(); i++) {
      bool_affected.push_back(false);
    }

    bk_state_time.resize(p_neu_pop->n_neurons());

    // `spike_events' holds spike events between `bk_state_time' and `t_end'
    TySpikeEventVec spike_events;
    // for fixed step, should contain no event
    TySpikeEventVec tmp_spike_events;

    // Save states at time `bk_state_time', then do a test step.
    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    for (double &tp : bk_state_time) tp = t;
    bk_neu_state = p_neu_pop->GetDymState();
    poisson_time_vec.SaveIdxAndClean();
    NextStepNoInteract(p_neu_pop, spike_events, t_end - t);

    // Evolve along the accurate order of spikes.
    // Each loop deals with one spike.
    while (!spike_events.empty()) {
      // Find out the first spike.
      auto heading_it = std::min_element(spike_events.begin(), spike_events.end());
      heading_spike_event = *heading_it;
      spike_events.erase(heading_it);     // "pop" the event
      dbg_printf("Dealing neuron %d spike at t = %f. Spikes left: %lu\n",
                 heading_spike_event.id, heading_spike_event.time, spike_events.size());

      // Find the index of affected neurons. Including the spiking one.
      for (const int &id : ids_affected)  // remove hash of last loop
        bool_affected[id] = false;
      ids_affected.clear();
      dbg_printf("  Affected neurons:");
      const auto &net = p_neu_pop->GetNeuronalParamsPtr()->net;
      for (SparseMat::InnerIterator it(net, heading_spike_event.id); it; ++it) {
        ids_affected.push_back(it.row());
        bool_affected[it.row()] = true;
        dbg_printf("  %d", it.row());
      }
      ids_affected.push_back(heading_spike_event.id);
      bool_affected[heading_spike_event.id] = true;
      dbg_printf(".\n");

      // Roll back the affected neurons.
      p_neu_pop->ScatterCopy(bk_neu_state, ids_affected);
      poisson_time_vec.RestoreIdx(ids_affected);
      spike_events.erase(
        remove_if(spike_events.begin(), spike_events.end(),
          [&](const TySpikeEvent &e) {
            return bool_affected[e.id] &&
              !(bk_state_time[e.id]==heading_spike_event.time
                && e.time==heading_spike_event.time);
          }),
        spike_events.end());

      // Evolve the affected neurons to the time of first spike. (fix step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(p_neu_pop, bk_state_time, ids_affected,
        tmp_spike_events, heading_spike_event.time);
      // deal with unexpected spikes
      bool b_heading_spike_pushed = false;
      for (size_t i = 0; i < tmp_spike_events.size(); i++) {
        if (tmp_spike_events[i].id == heading_spike_event.id) {
          b_heading_spike_pushed = true;
        } else {
          // possible extra spikes due to inaccurate computation
          cerr << "Unexpected spike (spike before \"first\" spike):  ["
            << i << "] id = " << tmp_spike_events[i].id
            << " time = " << tmp_spike_events[i].time << "\n";
          spike_events.emplace_back(heading_spike_event.time,
              tmp_spike_events[i].id);
        }
      }
      if (!b_heading_spike_pushed) {
        p_neu_pop->ForceReset(heading_spike_event.id);
      }

      // Perform synaptic interactions.
      p_neu_pop->SynapticInteraction(heading_spike_event);
      ras.emplace_back(heading_spike_event);
      vec_n_spike[heading_spike_event.id]++;

      // Save roll back base.
      for (const int &id : ids_affected)
        bk_state_time[id] = heading_spike_event.time;
      bk_neu_state.ScatterCopy(p_neu_pop->GetDymState(), ids_affected);
      poisson_time_vec.SaveIdxAndClean(ids_affected);

      // Evolve the affected neurons to t_end. (test step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(p_neu_pop, bk_state_time, ids_affected,
        spike_events, t_end);
    }
    t = t_end;
  }
};

class NeuronSimulatorExactSpikeOrderSparse2
: public NeuronSimulatorExactSpikeOrderSparse
{
public:

  NeuronSimulatorExactSpikeOrderSparse2(const TyNeuronalParams &_pm, double _dt)
    :NeuronSimulatorExactSpikeOrderSparse(_pm, _dt)
  {
  }

  // Evolve the whole system one dt, with interaction.
  // Save spike events in this `dt' in `ras',
  // and add spike count to `vec_n_spike'.
  MACRO_NO_INLINE void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;

    struct TySpikeEvent heading_spike_event(qNaN, -1);
    std::vector<int> ids_affected;             // index of affected neurons
    std::vector<bool> bool_affected;           // hash for `ids_affected'
    std::vector<bool> bool_fired;  // assume one neuron can not fire twice in one dt
    for (int i = 0; i < p_neu_pop->n_neurons(); i++) {
      bool_affected.push_back(false);
      bool_fired.push_back(false);
    }

    // `bk_neu_state' is state at time `bk_state_time'
    TyArrVals bk_state_time(p_neu_pop->n_neurons());
    struct TyNeuronalDymState bk_neu_state;

    // `spike_events' holds spike events between `bk_state_time' and `t_end'
    TySpikeEventVec spike_events;
    // for fixed step, should contain no event
    TySpikeEventVec tmp_spike_events;

    // Save states at time `bk_state_time', then do a test step.
    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    for (double &tp : bk_state_time) tp = t;
    bk_neu_state = p_neu_pop->GetDymState();
    poisson_time_vec.SaveIdxAndClean();
    NextStepNoInteract(p_neu_pop, spike_events, t_end - t);

    // Evolve in the accurate order of spikes.
    // Each loop deals with one spike.
    while (!spike_events.empty()) {
      // Find out the first spike.
      auto heading_it = std::min_element(spike_events.begin(), spike_events.end());
      heading_spike_event = *heading_it;
      spike_events.erase(heading_it);     // "pop" the event
      bool_fired[heading_spike_event.id] = true;
      dbg_printf("Dealing neuron %d spike at t = %f. Spikes left: %lu\n",
                 heading_spike_event.id, heading_spike_event.time, spike_events.size());

      // Find the index of affected neurons. Not including the spiking one.
      for (const int &id : ids_affected)  // remove hash of last loop
        bool_affected[id] = false;
      ids_affected.clear();
      dbg_printf("  Affected neurons:");
      const auto &net = p_neu_pop->GetNeuronalParamsPtr()->net;
      for (SparseMat::InnerIterator it(net, heading_spike_event.id); it; ++it) {
        ids_affected.push_back(it.row());
        bool_affected[it.row()] = true;
        dbg_printf("  %d", it.row());
      }
      dbg_printf(".\n");

      // Roll back the affected neurons.
      p_neu_pop->ScatterCopy(bk_neu_state, ids_affected);
      poisson_time_vec.RestoreIdx(ids_affected);
      spike_events.erase(
        remove_if(spike_events.begin(), spike_events.end(),
          [&](const TySpikeEvent &e) { return bool_affected[e.id]; }),
        spike_events.end());

      // Evolve the affected neurons to the time of first spike. (fix step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(p_neu_pop, bk_state_time, ids_affected,
        tmp_spike_events, heading_spike_event.time);
      // deal with unexpected spikes
      for (size_t i = 0; i < tmp_spike_events.size(); i++) {
        if (bool_fired[tmp_spike_events[i].id])
          continue;
        // possible extra spikes due to inaccurate computation
        cerr << "Unexpected spike (spike before \"first\" spike):  ["
          << i << "] id = " << tmp_spike_events[i].id
          << " time = " << tmp_spike_events[i].time << "\n";
        spike_events.emplace_back(heading_spike_event.time,
            tmp_spike_events[i].id);
      }

      // Perform synaptic interactions.
      p_neu_pop->SynapticInteraction(heading_spike_event);
      ras.emplace_back(heading_spike_event);
      vec_n_spike[heading_spike_event.id]++;

      // Save roll back base.
      for (const int &id : ids_affected)
        bk_state_time[id] = heading_spike_event.time;
      bk_neu_state.ScatterCopy(p_neu_pop->GetDymState(), ids_affected);
      poisson_time_vec.SaveIdxAndClean(ids_affected);

      // Evolve the affected neurons to t_end. (test step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(p_neu_pop, bk_state_time, ids_affected,
        tmp_spike_events, t_end);
      for (size_t i = 0; i < tmp_spike_events.size(); i++) {
        if (bool_fired[tmp_spike_events[i].id])
          continue;
        spike_events.emplace_back(tmp_spike_events[i]);
      }
    }
    t = t_end;
  }
};

#endif
