#ifndef HEADER_SIMULATOR_EXACT_ORDER
#define HEADER_SIMULATOR_EXACT_ORDER

#include "common_header.h"
#include "poisson_generator.h"
#include "neuron_system_utils.h"

/**
  Solver for pulsed coupled neuron model:
    After a test step (without synaptic interaction), evolve the system
    to the time of the first spike in this delta-t step, perform the synaptic
    interaction. Then start another test step for the remaining interval,
    recursively until there is no spike in the remaining interval.
  The template parameter TyNeuronModel fully describes the single neuron dynamics.
*/
template<typename TyNeuronModel>
class NeuronSimulatorExactSpikeSeqEvolve
{
public:
  TyNeuronModel neuron_model;
  struct TyNeuronalParams pm;
  struct TyNeuronalDymState<TyNeuronModel> neu_state;
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  NeuronSimulatorExactSpikeSeqEvolve(const TyNeuronModel &_neuron_model, const TyNeuronalParams &_pm, double _dt)
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
      double spike_time_local = std::numeric_limits<double>::quiet_NaN();
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
  // Time cost:
  //   # calls to NextStepSingleNeuronQuiet(): p*T/dt + p*pr*T + 2*p*p*fr*T
  //   synaptic interaction: (# of edges)*fr*T.
  // where the fr is mean firing rate over all neurons.
  __attribute__ ((noinline)) void NextDt(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;
    struct TyNeuronalDymState<TyNeuronModel> bk_neu_state;
    struct TySpikeEvent heading_spike_event(std::numeric_limits<double>::quiet_NaN(), -1);
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
        poisson_time_vec.RestoreIdx();             // replay the poisson events
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
            cerr << "Unexpected spike (spike before \"first\" spike):  [" << i << "] id = " << spike_events[i].id << " time = " << spike_events[i].time << "\n";
          }
          // Here does not do `ras.emplace_back(spike_events[i])' because
          // effectively the spike interaction are done at the same time.
          ras.emplace_back(heading_spike_event.time, spike_events[i].id);
          vec_n_spike[spike_events[i].id]++;
          SynapticInteraction(neu_state, spike_events[i]);
        }
        // Force the neuron to spike, if not already.
        if (!b_heading_spike_pushed) {
          neuron_model.VoltHandReset(neu_state.StatePtr(heading_spike_event.id));
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
      if (neu_state.dym_vals(j, neuron_model.id_V) > neuron_model.V_threshold) {
        neu_state.dym_vals(j, neuron_model.id_V) = 0;
        neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();
      }
    }
  }
};

#endif