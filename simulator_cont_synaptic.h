#ifndef HEADER_SIMULATOR_CONT_SYNAPTIC
#define HEADER_SIMULATOR_CONT_SYNAPTIC

#include "common_header.h"
#include "poisson_generator.h"
#include "neuron_system_utils.h"
#include "neuron_population.h"
#include "simulator_base.h"

class NeuronSimulatorCont :public NeuronSimulatorBase
{
public:
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  NeuronSimulatorCont(const TyNeuronalParams &pm, double _dt)
  {
    dt = _dt;
    t = 0;
    poisson_time_vec.Init(pm.arr_pr, t);
  }

  double GetT() const
  {
    return t;
  }

  TySpikeEventVec tmp_poisson_events;

public:
  void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector<size_t> &vec_n_spike) override
  {
    double t_local = t;
    double t_step_end = t + dt;
//    printf("From t = %.16e to %.16e\n", t, t_step_end);

    poisson_time_vec.SaveIdxAndClean();  // Clean Poisson queue

    // Get all Poisson events in this dt.
    tmp_poisson_events.clear();
    for (int j = 0; j < p_neu_pop->n_neurons(); j++) {
      TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
      double pr = p_neu_pop->GetNeuronalParamsPtr()->arr_pr[j];
      while (poisson_time_seq.Front() < t_step_end) {
        tmp_poisson_events.emplace_back(poisson_time_seq.Front(), j);
        poisson_time_seq.PopAndFill(pr);
      }
    }
    sort(tmp_poisson_events.begin(), tmp_poisson_events.end());

    size_t ras_old_size = ras.size();

    // Evoluate one dt, loop over Poisson events in this dt.
    for (size_t i = 0; i < tmp_poisson_events.size(); i++) {
//      printf("  P: to %d @ t = %.16e\n", tmp_poisson_events[i].id, tmp_poisson_events[i].time);
      p_neu_pop->NoInteractDt(tmp_poisson_events[i].time - t_local, t_local, ras);
      t_local = tmp_poisson_events[i].time;
      p_neu_pop->InjectPoissonE(tmp_poisson_events[i].id);
    }
    p_neu_pop->NoInteractDt(t_step_end - t_local, t_local, ras);

    for (auto it = ras.begin()+ras_old_size; it != ras.end(); it++) {
      vec_n_spike[it->id]++;
    }

    t += dt;
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

#endif
