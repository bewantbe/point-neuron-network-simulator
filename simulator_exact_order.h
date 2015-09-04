#ifndef HEADER_SIMULATOR_EXACT_ORDER
#define HEADER_SIMULATOR_EXACT_ORDER

#include "common_header.h"
#include "poisson_generator.h"
#include "neuron_system_utils.h"

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
template<typename TyNeuronModel>
class NeuronSimulatorExactSpikeOrder
{
public:
  TyNeuronModel neuron_model;
  struct TyNeuronalParams pm;
  struct TyNeuronalDymState<TyNeuronModel> neu_state;
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  NeuronSimulatorExactSpikeOrder(const TyNeuronModel &_neuron_model, const TyNeuronalParams &_pm, double _dt)
  :neuron_model(_neuron_model), pm(_pm), neu_state(_pm)
  {
    dt = _dt;
    t = 0;
    poisson_time_vec.Init(pm.arr_pr, t);
  }

protected:
  __attribute__ ((noinline)) void SynapticInteraction(struct TyNeuronalDymState<TyNeuronModel> &e_neu_state, const TySpikeEvent &se) const
  {
    if (se.id < pm.n_E) {
      // Excitatory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), neuron_model.id_gEInject) += it.value() * pm.scee;
        } else {
          e_neu_state.dym_vals(it.row(), neuron_model.id_gEInject) += it.value() * pm.scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), neuron_model.id_gIInject) += it.value() * pm.scei;
        } else {
          e_neu_state.dym_vals(it.row(), neuron_model.id_gIInject) += it.value() * pm.scii;
        }
      }
    }
  }

  // Evolve all neurons without synaptic interaction
  __attribute__ ((noinline)) void NextStepNoInteract(struct TyNeuronalDymState<TyNeuronModel> &tmp_neu_state,
              TySpikeEventVec &spike_events, double dt)
  {
    double t_step_end = t + dt;
    dbg_printf("----- Trying t = %f .. %f -------\n", t, t_step_end);
    for (int j = 0; j < pm.n_total(); j++) {
      //! tmp_neu_state.dym_vals must be Row major !
      TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
      double *dym_val = tmp_neu_state.StatePtr(j);
      double t_local = t;
      double spike_time_local = qNaN;
      while (poisson_time_seq.Front() < t_step_end) {
        dbg_printf("Neuron %d receive %lu-th (%lu) Poisson input at %f\n",
                   j, poisson_time_seq.id_seq, poisson_time_seq.size(),
                   poisson_time_seq.Front());
        dbg_printf("Time from %f to %f\n", t_local, poisson_time_seq.Front());
        neuron_model.NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
                            spike_time_local, poisson_time_seq.Front() - t_local);
        if (!std::isnan(spike_time_local)) {
          spike_events.emplace_back(t_local + spike_time_local, j);
        }
        t_local = poisson_time_seq.Front();
        dym_val[neuron_model.id_gEInject] += pm.arr_ps[j];  // Add Poisson input
        poisson_time_seq.PopAndFill(pm.arr_pr[j]);
      }
      dbg_printf("Time from %f to %f\n", t_local, t_step_end);
      neuron_model.NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
                          spike_time_local, t_step_end - t_local);
      if (!std::isnan(spike_time_local)) {
        spike_events.emplace_back(t_local + spike_time_local, j);
      }
      dbg_printf("Neuron %d @t=%f end state %f, %f, %f\n", j, t_step_end, dym_val[0], dym_val[1], dym_val[2]);
    }
  }

public:
  // Evolve the whole system one dt, with interaction
  __attribute__ ((noinline)) void NextDt(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;
    struct TyNeuronalDymState<TyNeuronModel> bk_neu_state;
    struct TySpikeEvent heading_spike_event(qNaN, -1);
    TySpikeEventVec spike_events;

    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    poisson_time_vec.SaveIdxAndClean();
    while (true) {
      spike_events.clear();
      bk_neu_state = neu_state;
      NextStepNoInteract(neu_state, spike_events, t_end - t);  // Try evolve till the step end
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
        neu_state = bk_neu_state;
        NextStepNoInteract(neu_state, spike_events, heading_spike_event.time - t);
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
          SynapticInteraction(neu_state, spike_events[i]);
        }
        // Force the neuron to spike, if not already.
        if (!b_heading_spike_pushed) {
          neuron_model.VoltHandReset(
              neu_state.StatePtr(heading_spike_event.id));
          neu_state.time_in_refractory[heading_spike_event.id] =
            std::numeric_limits<double>::min();
          ras.emplace_back(heading_spike_event);
          vec_n_spike[heading_spike_event.id]++;
          SynapticInteraction(neu_state, heading_spike_event);
        }
        t = heading_spike_event.time;
        poisson_time_vec.SaveIdxAndClean();
      }
    }
  }

  // Test and correct the neuron voltage (e.g. for imported data).
  void SaneTestVolt()
  {
    for (int j = 0; j < pm.n_total(); j++) {
      if (neu_state.dym_vals(j, neuron_model.id_V)
          > neuron_model.V_threshold) {
        neuron_model.VoltHandReset(neu_state.StatePtr(j));
        neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();
      }
    }
  }

  // Test if the calculation blow up
  void SaneTestState()
  {
    for (int j = 0; j < pm.n_total(); j++) {
      if ( !(fabs(neu_state.dym_vals(j, neuron_model.id_V))<500) ) {  // 500mV
        cerr << "\nNon-finite state value! V = " << neu_state.dym_vals(j, neuron_model.id_V)
          << "\n  Possible reason: time step too large.\n\n";
        throw "Non-finite state value!";
      }
    }
  }

  void PrintfState(const char *const rec_path, const struct TyNeuronalDymState<TyNeuronModel> &neu_state)
  {
    std::ofstream fout(rec_path, std::ios_base::app);
    fout << "\n";
    fout << "State in t = " << t << "\n";
    for (int j = 0; j < pm.n_total(); j++) {
      fout << "neu[" << std::setw(2) << j << "] = ";
      fout << setiosflags(std::ios::fixed) << std::setprecision(3);
      for (int k = 0; k < TyNeuronModel::n_var; k++) {
        fout << std::setw(7) << neu_state.dym_vals(j, k);
      }
      fout << "\n";
    }
    fout << std::endl;
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
template<typename TyNeuronModel>
class NeuronSimulatorExactSpikeOrderSparse
: public NeuronSimulatorExactSpikeOrder<TyNeuronModel>
{
  // Declare member variable/function for this template class.
  typedef NeuronSimulatorExactSpikeOrder<TyNeuronModel> NSE;
  using NSE::NextStepNoInteract;
  using NSE::SynapticInteraction;

public:
  using NSE::t;
  using NSE::dt;
  using NSE::pm;
  using NSE::neu_state;
  using NSE::neuron_model;
  using NSE::poisson_time_vec;

  NeuronSimulatorExactSpikeOrderSparse(const TyNeuronModel &_neuron_model,
    const TyNeuronalParams &_pm, double _dt)
    :NeuronSimulatorExactSpikeOrder<TyNeuronModel>(_neuron_model, _pm, _dt)
  {
  }

protected:
  // Evolve all neurons without synaptic interaction
  __attribute__ ((noinline)) void NextStepNoInteractToTime(
      struct TyNeuronalDymState<TyNeuronModel> &tmp_neu_state,
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
      double *dym_val = tmp_neu_state.StatePtr(j);
      double t_local = bk_state_time[j];
      double spike_time_local = qNaN;
      while (poisson_time_seq.Front() < t_step_end) {
        dbg_printf("  receive %lu-th (%lu) Poisson input at %f\n",
                   poisson_time_seq.id_seq, poisson_time_seq.size(),
                   poisson_time_seq.Front());
        dbg_printf("  time from %f to %f\n", t_local, poisson_time_seq.Front());
        neuron_model.NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
                            spike_time_local, poisson_time_seq.Front() - t_local);
        if (!std::isnan(spike_time_local)) {
          spike_events.emplace_back(t_local + spike_time_local, j);
        }
        t_local = poisson_time_seq.Front();
        dym_val[neuron_model.id_gEInject] += pm.arr_ps[j];  // Add Poisson input
        poisson_time_seq.PopAndFill(pm.arr_pr[j]);
      }
      dbg_printf("  time from %f to %f\n", t_local, t_step_end);
      neuron_model.NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
                          spike_time_local, t_step_end - t_local);
      if (!std::isnan(spike_time_local)) {
        spike_events.emplace_back(t_local + spike_time_local, j);
      }
      dbg_printf("  t=%f end state: %f, %f, %f\n", t_step_end, dym_val[0], dym_val[1], dym_val[2]);
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
  struct TyNeuronalDymState<TyNeuronModel> bk_neu_state;

  __attribute__ ((noinline)) void NextDt(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;

    struct TySpikeEvent heading_spike_event(qNaN, -1);
    ids_affected.clear();
    bool_affected.clear();
    for (int i = 0; i < pm.n_total(); i++) {
      bool_affected.push_back(false);
    }

    bk_state_time.resize(pm.n_total());

    // `spike_events' holds spike events between `bk_state_time' and `t_end'
    TySpikeEventVec spike_events;
    // for fixed step, should contain no event
    TySpikeEventVec tmp_spike_events;

    // Save states at time `bk_state_time', then do a test step.
    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    for (double &tp : bk_state_time) tp = t;
    bk_neu_state = neu_state;
    poisson_time_vec.SaveIdxAndClean();
    NextStepNoInteract(neu_state, spike_events, t_end - t);

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
      for (SparseMat::InnerIterator it(pm.net, heading_spike_event.id); it; ++it) {
        ids_affected.push_back(it.row());
        bool_affected[it.row()] = true;
        dbg_printf("  %d", it.row());
      }
      ids_affected.push_back(heading_spike_event.id);
      bool_affected[heading_spike_event.id] = true;
      dbg_printf(".\n");

      // Roll back the affected neurons.
      neu_state.ScatterCopy(bk_neu_state, ids_affected);
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
      NextStepNoInteractToTime(neu_state, bk_state_time, ids_affected,
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
        neuron_model.VoltHandReset(
            neu_state.StatePtr(heading_spike_event.id));
        neu_state.time_in_refractory[heading_spike_event.id] =
          std::numeric_limits<double>::min();
      }

      // Perform synaptic interactions.
      SynapticInteraction(neu_state, heading_spike_event);
      ras.emplace_back(heading_spike_event);
      vec_n_spike[heading_spike_event.id]++;

      // Save roll back base.
      for (const int &id : ids_affected)
        bk_state_time[id] = heading_spike_event.time;
      bk_neu_state.ScatterCopy(neu_state, ids_affected);
      poisson_time_vec.SaveIdxAndClean(ids_affected);

      // Evolve the affected neurons to t_end. (test step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(neu_state, bk_state_time, ids_affected,
        spike_events, t_end);
    }
    t = t_end;
  }
};

template<typename TyNeuronModel>
class NeuronSimulatorExactSpikeOrderSparse2
: public NeuronSimulatorExactSpikeOrderSparse<TyNeuronModel>
{
  // Declare member variable/function for this template class.
  typedef NeuronSimulatorExactSpikeOrderSparse<TyNeuronModel> NSE;
  using NSE::NextStepNoInteract;
  using NSE::SynapticInteraction;
  using NSE::NextStepNoInteractToTime;

public:
  using NSE::t;
  using NSE::dt;
  using NSE::pm;
  using NSE::neu_state;
  using NSE::neuron_model;
  using NSE::poisson_time_vec;

  NeuronSimulatorExactSpikeOrderSparse2(const TyNeuronModel &_neuron_model,
    const TyNeuronalParams &_pm, double _dt)
    :NeuronSimulatorExactSpikeOrderSparse<TyNeuronModel>(_neuron_model, _pm, _dt)
  {
  }

  // Evolve the whole system one dt, with interaction.
  // Save spike events in this `dt' in `ras',
  // and add spike count to `vec_n_spike'.
  __attribute__ ((noinline)) void NextDt(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;

    struct TySpikeEvent heading_spike_event(qNaN, -1);
    std::vector<int> ids_affected;             // index of affected neurons
    std::vector<bool> bool_affected;           // hash for `ids_affected'
    std::vector<bool> bool_fired;  // assume one neuron can not fire twice in one dt
    for (int i = 0; i < pm.n_total(); i++) {
      bool_affected.push_back(false);
      bool_fired.push_back(false);
    }

    // `bk_neu_state' is state at time `bk_state_time'
    TyArrVals bk_state_time(pm.n_total());
    struct TyNeuronalDymState<TyNeuronModel> bk_neu_state;

    // `spike_events' holds spike events between `bk_state_time' and `t_end'
    TySpikeEventVec spike_events;
    // for fixed step, should contain no event
    TySpikeEventVec tmp_spike_events;

    // Save states at time `bk_state_time', then do a test step.
    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    for (double &tp : bk_state_time) tp = t;
    bk_neu_state = neu_state;
    poisson_time_vec.SaveIdxAndClean();
    NextStepNoInteract(neu_state, spike_events, t_end - t);

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
      for (SparseMat::InnerIterator it(pm.net, heading_spike_event.id); it; ++it) {
        ids_affected.push_back(it.row());
        bool_affected[it.row()] = true;
        dbg_printf("  %d", it.row());
      }
      dbg_printf(".\n");

      // Roll back the affected neurons.
      neu_state.ScatterCopy(bk_neu_state, ids_affected);
      poisson_time_vec.RestoreIdx(ids_affected);
      spike_events.erase(
        remove_if(spike_events.begin(), spike_events.end(),
          [&](const TySpikeEvent &e) { return bool_affected[e.id]; }),
        spike_events.end());

      // Evolve the affected neurons to the time of first spike. (fix step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(neu_state, bk_state_time, ids_affected,
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
      SynapticInteraction(neu_state, heading_spike_event);
      ras.emplace_back(heading_spike_event);
      vec_n_spike[heading_spike_event.id]++;

      // Save roll back base.
      for (const int &id : ids_affected)
        bk_state_time[id] = heading_spike_event.time;
      bk_neu_state.ScatterCopy(neu_state, ids_affected);
      poisson_time_vec.SaveIdxAndClean(ids_affected);

      // Evolve the affected neurons to t_end. (test step)
      tmp_spike_events.clear();
      NextStepNoInteractToTime(neu_state, bk_state_time, ids_affected,
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
