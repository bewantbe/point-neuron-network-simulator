#ifndef HEADER_SIMULATOR_BASE
#define HEADER_SIMULATOR_BASE

#include "common_header.h"
#include "neuron_system_utils.h"

class NeuronSimulatorBase
{
public:
  virtual void NextDt(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike) = 0;
  virtual void SaneTestVolt() = 0;
  virtual void SaneTestState() = 0;
  virtual const TyNeuronalDymState & GetNeuState() const = 0;
  virtual TyNeuronalDymState & GetNeuState() = 0;
  virtual const TyPoissonTimeVec & Get_poisson_time_vec() const = 0;
  virtual TyPoissonTimeVec & Get_poisson_time_vec() = 0;
};

#endif
