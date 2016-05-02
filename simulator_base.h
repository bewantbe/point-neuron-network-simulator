#ifndef HEADER_SIMULATOR_BASE
#define HEADER_SIMULATOR_BASE

#include "common_header.h"
#include "poisson_generator.h"
#include "neuron_system_utils.h"
#include "neuron_population.h"

class NeuronSimulatorBase
{
public:
  virtual void NextDt(NeuronPopulationBase * p_neu_pop,
      TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike) = 0;
  virtual void SaneTestVolt() = 0;
  virtual void SaneTestState() = 0;
  virtual const TyPoissonTimeVec & Get_poisson_time_vec() const = 0;
  virtual TyPoissonTimeVec & Get_poisson_time_vec() = 0;
};

#endif
