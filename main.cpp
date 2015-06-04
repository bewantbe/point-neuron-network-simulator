// g++ -g -O0 -std=c++11 -Wall main.cpp math_helper.cpp -lboost_program_options -o bin/vec_IFsimu
//Could faster:
// g++ -g -O2 -falign-functions=16 -falign-loops=16

#define NDEBUG

#include <cassert>
#include "common_header.h"
#include "math_helper.h"

// Fast function calculation for exp(x)
#include "fmath.hpp"
#define exp(x) fmath::expd(x)

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//const bool g_b_debug = false;
#ifdef NDEBUG
#define dbg_printf(...) ((void)0);
#else
#define dbg_printf printf
#endif

using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMat;
typedef std::vector<double> TyArrVals;

// list of neuron models
enum eNeuronModel
{
  LIF_G,
  LIF_GH,
  HH
};

struct Ty_LIF_G
{
  // The neuron model named LIF-G in this code.
  // This is the Leaky Integrate-and-Fire model with jump conductance.
  double Vot_Threshold   = 1.0;  // voltages are in dimensionless unit
  double Vot_Reset       = 0.0;
  double Vot_Leakage     = 0.0;
  double Vot_Excitatory  = 14.0/3.0;
  double Vot_Inhibitory  = -2.0/3.0;
  double Con_Leakage     = 0.05;  // ms^-1
  double Time_ExCon      = 2.0;   // ms
  double Time_InCon      = 5.0;   // ms
  double Time_Refractory = 2.0;   // ms
  static const int n_var = 3;  // number of dynamical variables
  static const int id_V  = 0;  // index of V variable
  static const int id_gE = 1;  // index of gE variable
  static const int id_gI = 2;  // index of gI variable
  static const int id_gEInject = id_gE;  // index of gE injection variable
  static const int id_gIInject = id_gI;  // index of gI injection variable

  // Evolve conductance only
  void NextDtConductance(double *dym_val, double dt) const
  {
    dym_val[id_gE] *= exp(-dt / Time_ExCon);
    dym_val[id_gI] *= exp(-dt / Time_InCon);
  }

  // Get instantaneous dv/dt for current dynamical state
  double GetDv(const double *dym_val) const
  {
    return - Con_Leakage * (dym_val[id_V] - Vot_Leakage)
           - dym_val[id_gE] * (dym_val[id_V] - Vot_Excitatory)
           - dym_val[id_gI] * (dym_val[id_V] - Vot_Inhibitory);
  }

  // Evolve the state `dym_val' a `dt' forward,
  // using classical Runge–Kutta 4-th order scheme for voltage.
  // Conductance will evolve using the exact formula.
  // Return derivative k1 at t_n, for later interpolation.
  __attribute__ ((noinline)) double DymInplaceRK4(double *dym_val, double dt) const
  {

    double v_n = dym_val[id_V];
    double k1, k2, k3, k4;
    double exp_E = exp(-0.5 * dt / Time_ExCon);
    double exp_I = exp(-0.5 * dt / Time_InCon);

    // k1 = f(t_n, y_n)
    k1 = GetDv(dym_val);

    // y_n + 0.5*k1*dt
    dym_val[id_V ] = v_n + 0.5 * dt * k1;
    dym_val[id_gE] *= exp_E;
    dym_val[id_gI] *= exp_I;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k2 = GetDv(dym_val);

    // y_n + 0.5*k2*dt
    dym_val[id_V ] = v_n + 0.5 * dt * k2;
    // k3 = f(t+dt/2, y_n + 0.5*k2*dt)
    k3 = GetDv(dym_val);

    // y_n + k3*dt
    dym_val[id_V ] = v_n + dt * k3;
    dym_val[id_gE] *= exp_E;
    dym_val[id_gI] *= exp_I;
    // k4 = f(t+dt, y_n + k3*dt)
    k4 = GetDv(dym_val);

    dym_val[id_V ] = v_n + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);

    return k1;
  }
};

struct Ty_LIF_GH
{
  // The neuron model named LIF-GH in this code.
  // This is the Leaky Integrate-and-Fire model with order 1 smoothed conductance.
  double Vot_Threshold   = 1.0;  // voltages are in dimensionless unit
  double Vot_Reset       = 0.0;
  double Vot_Leakage     = 0.0;
  double Vot_Excitatory  = 14.0/3.0;
  double Vot_Inhibitory  = -2.0/3.0;
  double Con_Leakage     = 0.05;  // ms^-1
  double Time_ExCon      = 2.0;   // ms
  double Time_ExConR     = 0.5;   // ms
  double Time_InCon      = 5.0;   // ms
  double Time_InConR     = 0.8;   // ms
  double Time_Refractory = 2.0;   // ms
  static const int n_var    = 5;  // number of dynamical variables
  static const int id_V     = 0;  // index of V variable
  static const int id_gE    = 1;  // index of gE variable
  static const int id_gI    = 2;  // index of gI variable
  static const int id_gE_s1 = 3;  // index of derivative of gE
  static const int id_gI_s1 = 4;  // index of derivative of gI
  static const int id_gEInject = id_gE_s1;  // index of gE injection variable
  static const int id_gIInject = id_gI_s1;  // index of gI injection variable

  // Evolve conductance only
  void NextDtConductance(double *dym_val, double dt) const
  {
    /**
      ODE:
        g '[t] == -1/tC  * g [t] + gR[t]
        gR'[t] == -1/tCR * gR[t]
      Solution:
        g [t] = exp(-t/tC) * g[0] + (exp(-t/tC) - exp(-t/tCR)) * gR[0] * tC * tCR / (tC - tCR)
        gR[t] = exp(-t/tCR) * gR[0]
      Another form (hopefully more accurate, note exp(x)-1 is expm1(x) ):
        g [t] = exp(-t/tC) * (g[0] + gR[0] * (exp((1/tC-1/tCR)*t) - 1) / (1/tC - 1/tCR))
      Or
        g [t] = exp(-t/tC) * g[0] + exp(-t/tCR) * gR[0] * (exp((1/tCR-1/tC)*t) - 1) / (1/tCR-1/tC)
    */
    // Excitatory
    double expC  = exp(-dt / Time_ExCon);
    double expCR = exp(-dt / Time_ExConR);
    dym_val[id_gE] = expC*dym_val[id_gE] + (expC - expCR) * dym_val[id_gE_s1] * Time_ExCon * Time_ExConR / (Time_ExCon - Time_ExConR);
    dym_val[id_gE_s1] *= expCR;
    // Inhibitory
    expC  = exp(-dt / Time_InCon);
    expCR = exp(-dt / Time_InConR);
    dym_val[id_gI] = expC*dym_val[id_gI] + (expC - expCR) * dym_val[id_gI_s1] * Time_InCon * Time_InConR / (Time_InCon - Time_InConR);
    dym_val[id_gI_s1] *= expCR;
  }

  // Get instantaneous dv/dt for current dynamical state
  double GetDv(const double *dym_val) const
  {
    return - Con_Leakage * (dym_val[id_V] - Vot_Leakage)
           - dym_val[id_gE] * (dym_val[id_V] - Vot_Excitatory)
           - dym_val[id_gI] * (dym_val[id_V] - Vot_Inhibitory);
  }

  // Evolve the state `dym_val' a `dt' forward,
  // using classical Runge–Kutta 4-th order scheme for voltage.
  // Conductance will evolve using the exact formula.
  // Return derivative k1 at t_n, for later interpolation.
  __attribute__ ((noinline)) double DymInplaceRK4(double *dym_val, double dt) const
  {
    double v_n = dym_val[id_V];
    double k1, k2, k3, k4;
    double expEC  = exp(-0.5 * dt / Time_ExCon);  // TODO: maybe cache this value?
    double expECR = exp(-0.5 * dt / Time_ExConR);
    double expIC  = exp(-0.5 * dt / Time_InCon);
    double expICR = exp(-0.5 * dt / Time_InConR);
    double gE_s_coef = (expEC - expECR) * Time_ExCon * Time_ExConR / (Time_ExCon - Time_ExConR);
    double gI_s_coef = (expIC - expICR) * Time_InCon * Time_InConR / (Time_InCon - Time_InConR);

    // k1 = f(t_n, y_n)
    k1 = GetDv(dym_val);

    // y_n + 0.5*k1*dt
    dym_val[id_V ] = v_n + 0.5 * dt * k1;
    dym_val[id_gE] = expEC * dym_val[id_gE] + gE_s_coef * dym_val[id_gE_s1];
    dym_val[id_gI] = expIC * dym_val[id_gI] + gI_s_coef * dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k2 = GetDv(dym_val);

    // y_n + 0.5*k2*dt
    dym_val[id_V ] = v_n + 0.5 * dt * k2;
    // k3 = f(t+dt/2, y_n + 0.5*k2*dt)
    k3 = GetDv(dym_val);

    // y_n + k3*dt
    dym_val[id_V ] = v_n + dt * k3;
    dym_val[id_gE] = expEC * dym_val[id_gE] + gE_s_coef * dym_val[id_gE_s1];
    dym_val[id_gI] = expIC * dym_val[id_gI] + gI_s_coef * dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    // k4 = f(t+dt, y_n + k3*dt)
    k4 = GetDv(dym_val);

    dym_val[id_V ] = v_n + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);

    return k1;
  }
};

// Model of HH neuron, with two DE for G
struct Ty_HH_GH
{
  double V_Na = 115;             // mV
  double V_K  = -12;
  double V_L  =  10.5;
  double G_Na = 120;             // mS cm^-2
  double G_K  =  36;
  double G_L  =   0.3;
  double V_gE =  65;             // mV
  double V_gI = -15;             // mV        
  double time_gE    = 0.5;       // ms
  double time_gE_s1 = 3.0;       // ms
  double time_gI    = 0.5;       // ms
  double time_gI_s1 = 7.0;       // ms
  double V_threshold = 15;       // mV, used to determine the spike timing
  static const int n_var = 8;
  static const int n_var_soma = 4;
  static const int id_V = 0;
  static const int id_m = 1;
  static const int id_h = 2;
  static const int id_n = 3;
  static const int id_gE = 4;
  static const int id_gI = 5;
  static const int id_gE_s1 = 6;
  static const int id_gI_s1 = 7;
  static const int id_gEInject = id_gE_s1;
  static const int id_gIInject = id_gI_s1;

  // should be inlined
  double GetDv(const double *dym_val) const
  {
    return
      -(dym_val[id_V]-V_Na) * G_Na
        * dym_val[id_h] * dym_val[id_m] * dym_val[id_m] * dym_val[id_m]
      -(dym_val[id_V]-V_K ) * G_K
        * dym_val[id_n] * dym_val[id_n] * dym_val[id_n] * dym_val[id_n]
      -(dym_val[id_V]-V_L ) * G_L
      -(dym_val[id_V]-V_gE) * dym_val[id_gE]
      -(dym_val[id_V]-V_gI) * dym_val[id_gI];

  }

  void ODE_RHS(const double *dym_val, double *dym_d_val) const
  {
    dym_d_val[id_V] = GetDv(dym_val);
    dym_d_val[id_m] = (1-dym_val[id_m])
        * ((0.1-0.01*dym_val[id_V])/(exp(1-0.1*dym_val[id_V])-1))
      - dym_val[id_m] * (0.125*exp(-dym_val[id_V]/80));
    dym_d_val[id_h] = (1-dym_val[id_h])
        * ((2.5-0.1*dym_val[id_V])/(exp(2.5-0.1*dym_val[id_V])-1))
      - dym_val[id_h] * (4*exp(-dym_val[id_V]/18));
    dym_d_val[id_n] = (1-dym_val[id_n]) * (0.07*exp(-dym_val[id_V]/20))
      - dym_val[id_n] / (exp(3-0.1*dym_val[id_V])+1);
    // No RHS for G
  }

  double DymInplaceRK4(double *dym_val, double dt) const
  {
    double expEC  = exp(-0.5 * dt / time_gE);
    double expECR = exp(-0.5 * dt / time_gE_s1);
    double expIC  = exp(-0.5 * dt / time_gI);
    double expICR = exp(-0.5 * dt / time_gI_s1);
    double gE_s_coef = (expEC - expECR) * time_gE * time_gE_s1 / (time_gE - time_gE_s1);
    double gI_s_coef = (expIC - expICR) * time_gI * time_gI_s1 / (time_gI - time_gE_s1);

    double k1[n_var_soma], k2[n_var_soma], k3[n_var_soma], k4[n_var_soma];
    double dym_val0[n_var_soma];
    Eigen::Map< Eigen::RowVectorXd > dym_val0_v(dym_val0, n_var_soma);
    Eigen::Map< Eigen::RowVectorXd > dym_val_v (dym_val , n_var_soma);
    Eigen::Map< Eigen::RowVectorXd > k1_v(k1, n_var_soma);
    Eigen::Map< Eigen::RowVectorXd > k2_v(k2, n_var_soma);
    Eigen::Map< Eigen::RowVectorXd > k3_v(k3, n_var_soma);
    Eigen::Map< Eigen::RowVectorXd > k4_v(k4, n_var_soma);

    dym_val0_v = dym_val_v;
    ODE_RHS(dym_val, k1);

    dym_val_v = dym_val0_v + 0.5*dt*k1_v;
    dym_val[id_gE] = expEC*dym_val[id_gE] + gE_s_coef*dym_val[id_gE_s1];
    dym_val[id_gI] = expIC*dym_val[id_gI] + gI_s_coef*dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    ODE_RHS(dym_val, k2);

    dym_val_v = dym_val0_v + 0.5*dt*k2_v;
    ODE_RHS(dym_val, k3);

    dym_val_v = dym_val0_v + dt*k3_v;
    dym_val[id_gE] = expEC*dym_val[id_gE] + gE_s_coef*dym_val[id_gE_s1];
    dym_val[id_gI] = expIC*dym_val[id_gI] + gI_s_coef*dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    ODE_RHS(dym_val, k4);

    dym_val_v = dym_val0_v + dt/6.0 * (k1_v + 2*k2_v + 2*k3_v + k4_v);

    return k1[id_V];
  };

  void NextStepSingleNeuronQuiet(double *dym_val, double &t_in_refractory,
    double &spike_time_local, double dt_local) const
  {
    t_in_refractory = 0;  // pretend no refractory period
    spike_time_local = std::numeric_limits<double>::quiet_NaN();
    double v0 = dym_val[id_V];
    double k1 = DymInplaceRK4(dym_val, dt_local);
    double &v1 = dym_val[id_V];
    if (v0 <= V_threshold && v1 > V_threshold) {
      spike_time_local = cubic_hermit_real_root(dt_local,
        v0, v1, k1, GetDv(dym_val), V_threshold);
    }
  }
};

// Parameters about neuronal system
struct TyNeuronalParams
{
  int n_E, n_I;    // Number of neurons, Excitatory and Inhibitory type
  SparseMat net;   // Adjacency matrix, net(i,j) means i is affect by j
  double scee, scie, scei, scii;
  TyArrVals arr_pr;   // Poisson input rate for each neuron
  TyArrVals arr_ps;   // Poisson input strength for each neuron
  // TODO: maybe add Inhibitory Poisson input?
  //TyArrVals arr_psi;

  int n_total() const
  { return n_E + n_I; }

  // should use this this function to initialize neurons
  void SetNumberOfNeurons(int _n_E, int _n_I)
  {
    n_E = _n_E;
    n_I = _n_I;
    arr_pr.resize(n_total());
    arr_ps.resize(n_total());
    net.resize(n_total(), n_total());
  }

  TyNeuronalParams(int _n_E, int _n_I)
  :scee(0), scie(0), scei(0), scii(0)
  {
    SetNumberOfNeurons(_n_E, _n_I);
  }
};

// Dynamical variables (V and G etc) for the neuron system
template<typename TyNeuronModel>
struct TyNeuronalDymState
{
  // Current dynamical states of neurons
  // Use RowMajor, so state variables of each neuron are continuous on memory
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dym_vals;
  // Time resided in refractory period for each neuron
  TyArrVals time_in_refractory;

  void Zeros()
  {
    dym_vals.setZero();
    for (auto &i : time_in_refractory) i = 0;
  }

  // Pointer to the state of j-th neuron
  double *StatePtr(int j)
  {
    return dym_vals.data() + j * dym_vals.cols();
  }

  TyNeuronalDymState(const TyNeuronalParams &pm)
  {
    dym_vals.resize(pm.n_total(), TyNeuronModel::n_var);
    time_in_refractory.resize(pm.n_total());
    Zeros();
  }

  TyNeuronalDymState()
  { }
};

//std::random_device rd;
std::mt19937 rand_eng(1);  // rand_eng(rd())

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

struct TySpikeEvent
{
  double time;
  int id;

  TySpikeEvent()
  {
    time = std::numeric_limits<double>::quiet_NaN();
    id = -1;
  }

  TySpikeEvent(double _time, int _id)
  : time(_time), id(_id)
  {}

  bool operator < (const TySpikeEvent &b) const
  { return time < b.time; }

  bool operator > (const TySpikeEvent &b) const
  { return time > b.time; }

  bool operator == (const TySpikeEvent &b) const
  { return time == b.time && id == b.id; }
};

typedef std::priority_queue<
    TySpikeEvent,
    std::vector<TySpikeEvent>,
    std::greater<TySpikeEvent> > TySpikeEventQueue;

typedef std::vector< TySpikeEvent > TySpikeEventVec;

typedef std::vector<double> TyInternalVec;

// Used to generate Poisson events for one neuron
class TyPoissonTimeSeq: public TyInternalVec
{
public:
  typedef size_t TyIdx;
  TyIdx id_seq = 0;      // point to current event

  TyPoissonTimeSeq() = default;
  TyPoissonTimeSeq(const TyPoissonTimeSeq &pts)
  : TyInternalVec(pts),
    id_seq(pts.id_seq)
  {}

  // Requirements: rate >= 0, t_until <= inf.
  void AddEventsUntilTime(double rate, double t_until)
  {
    assert(rate >= 0);
    std::exponential_distribution<> exp_dis(rate);
    while (back() < t_until) {
      push_back( back() + exp_dis(rand_eng) );
    }
  }

  void Init(double rate, double t0)
  {
    assert(rate >= 0);
    clear();
    std::exponential_distribution<> exp_dis(rate);
    push_back(t0 + exp_dis(rand_eng));
    id_seq = 0;
  }

  double Front() const
  {
    return operator[](id_seq);
  }

  void PopAndFill(double rate, bool auto_shrink = false)
  {
    assert(id_seq < size());
    id_seq++;
    if (id_seq == size()) {
      assert(size() > 0);
      if (std::isfinite(back())) {
        if (auto_shrink) {
          Init(rate, back());
        }
        AddEventsUntilTime(rate, back() + 12.0 / rate);
      } else {
        id_seq--;
      }
    }
  }

  void Shrink()
  {
    assert(id_seq < size());
    erase(begin(), begin()+id_seq);
    id_seq = 0;
  }
};

// Used to generate Poisson events for each neurons
class TyPoissonTimeVec: public std::vector<TyPoissonTimeSeq>
{
  std::vector<TyPoissonTimeSeq::TyIdx> id_seq_vec;  // point to current event
public:
  void Init(const TyArrVals &rate_vec, double t0)
  {
    // each TyPoissonTimeSeq should have at least one event
    resize(rate_vec.size());
    for (size_t j = 0; j < rate_vec.size(); j++) {
      operator[](j).Init(rate_vec[j], t0);
    }
    id_seq_vec.resize(rate_vec.size());
  }

  TyPoissonTimeVec(const TyArrVals &rate_vec, double t0 = 0)
  {
    Init(rate_vec, t0);
  }

  TyPoissonTimeVec() = default;  // enable all other default constructors

  void RestoreIdx()
  {
    for (size_t j = 0; j < size(); j++) {
      operator[](j).id_seq = id_seq_vec[j];
    }
  }

  void SaveIdxAndClean()
  {
    int j = 0;
    for (iterator it = begin(); it != end(); ++it, ++j) {
      if (it->size() - it->id_seq < it->id_seq / 7) {  // The factor here is non-critical
        it->Shrink();
        dbg_printf("  ( %d-th TyPoissonTimeSeq size shrinked )\n", j);
      }
      id_seq_vec[j] = it->id_seq;
    }
  }
};

/**
  Solver for pulsed coupled neuron model:
    After a test step (without synaptic interaction), then evolve the system
    to the time of the first spike in this delta-t step, perform the synaptic
    interaction. Then start another test step for the remaining interval,
    recursively until there is no spike in the remaining interval.
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

  // Evolve the ODE and note down the spike time, assuming no reset and no external input.
  // `spike_time_local' should be guaranteed to be with in [0, dt] or NAN.
  __attribute__ ((noinline)) void NextStepSingleNeuronContinuous(double *dym_val, double &spike_time_local, double dt) const
  {
    double v_n = dym_val[neuron_model.id_V];
    double k1  = neuron_model.DymInplaceRK4(dym_val, dt);

    if (v_n <= neuron_model.Vot_Threshold
        && dym_val[neuron_model.id_V ] > neuron_model.Vot_Threshold) {
      spike_time_local = cubic_hermit_real_root(dt,
        v_n, dym_val[neuron_model.id_V ],
        k1, neuron_model.GetDv(dym_val), neuron_model.Vot_Threshold);
    } else {
      if (v_n > 0.996 && k1>0) { // the v_n > 0.996 is for dt=0.5 ms, LIF,G model
        // Try capture some missing spikes that the intermediate value passes
        // threshold, but both ends are lower than threshold.
        // Get quadratic curve from value of both ends and derivative from left end
        // Return the maximum point as `t_max_guess'
        double c = v_n;
        double b = k1;
        double a = (dym_val[neuron_model.id_V ] - c - b*dt)/(dt*dt);
        double t_max_guess = -b/(2*a);
        // In LIF-G, it can guarantee that a<0 (concave),
        // hence t_max_guess > 0. But in LIF-G model, we still need to
        // check 0 < t_max_guess
        if (0 < t_max_guess && t_max_guess < dt
            && (b*b)/(-4*a)+c >= neuron_model.Vot_Threshold) {
          //dbg_printf("Rare event: mid-dt spike captured, guess time: %f\n", t_max_guess);
          dbg_printf("NextStepSingleNeuronContinuous(): possible mid-dt spike detected:\n");
          dbg_printf("  Guessed max time: %f, t = %f, dt = %f\n", t_max_guess, t, dt);
          // root should in [0, t_max_guess]
          spike_time_local = cubic_hermit_real_root(dt,
            v_n, dym_val[neuron_model.id_V ],
            k1, neuron_model.GetDv(dym_val), neuron_model.Vot_Threshold);
        } else {
          spike_time_local = std::numeric_limits<double>::quiet_NaN();
        }
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // Evolve single neuron as if no external input.
  // Return first spike time in `spike_time_local', if any.
  __attribute__ ((noinline)) void NextStepSingleNeuronQuiet(double *dym_val, double &t_in_refractory,
                           double &spike_time_local, double dt_local)
  {
    //! at most one spike allowed during this dt_local
    if (t_in_refractory == 0) {
      dbg_printf("NextStepSingleNeuronQuiet(): dt_local = %f:\n", dt_local);
      dbg_printf("NextStepSingleNeuronQuiet(): begin state=%f,%f,%f\n",
                 dym_val[0], dym_val[1], dym_val[2]);
      NextStepSingleNeuronContinuous(dym_val, spike_time_local, dt_local);
      dbg_printf("NextStepSingleNeuronQuiet(): end   state=%f,%f,%f\n",
                 dym_val[0], dym_val[1], dym_val[2]);
      if (!std::isnan(spike_time_local)) {
        // Add `numeric_limits<double>::min()' to make sure t_in_refractory > 0.
        t_in_refractory = dt_local - spike_time_local
                          + std::numeric_limits<double>::min();
        dym_val[neuron_model.id_V] = neuron_model.Vot_Reset;
        dbg_printf("Reach threshold detected\n");
        if (t_in_refractory >= neuron_model.Time_Refractory) {
          // Short refractory period (< dt_local), neuron will active again.
          dt_local = t_in_refractory - neuron_model.Time_Refractory;
          t_in_refractory = 0;
          // Back to the activation time.
          neuron_model.NextDtConductance(dym_val, -dt_local);
          double spike_time_local_tmp;
          NextStepSingleNeuronContinuous(dym_val, spike_time_local_tmp, dt_local);
          if (!std::isnan(spike_time_local_tmp)) {
            cerr << "NextStepSingleNeuronQuiet(): Multiple spikes in one dt. Interaction dropped." << endl;
            cerr << "  dt_local = " << dt_local << '\n';
            cerr << "  spike_time_local = " << spike_time_local << '\n';
            cerr << "  t_in_refractory = " << t_in_refractory << '\n';
            cerr << "  dym_val[neuron_model.id_V] = " << dym_val[neuron_model.id_V] << '\n';
            cerr << "  dym_val[neuron_model.id_gE] = " << dym_val[neuron_model.id_gE] << '\n';
            throw "Multiple spikes in one dt.";
          }
        }
      }
    } else {
      // Neuron in refractory period.
      double dt_refractory_remain = neuron_model.Time_Refractory
                                 - t_in_refractory;
      if (dt_refractory_remain < dt_local) {
        // neuron will awake after dt_refractory_remain which is in this dt_local
        neuron_model.NextDtConductance(dym_val, dt_refractory_remain);
        assert( dym_val[neuron_model.id_V] == neuron_model.Vot_Reset );
        t_in_refractory = 0;
        NextStepSingleNeuronQuiet(dym_val, t_in_refractory,
                            spike_time_local, dt_local - dt_refractory_remain);
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
        neuron_model.NextDtConductance(dym_val, dt_local);
        t_in_refractory += dt_local;
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
        NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
                            spike_time_local, poisson_time_seq.Front() - t_local);
        if (!std::isnan(spike_time_local)) {
          spike_events.emplace_back(t_local + spike_time_local, j);
        }
        t_local = poisson_time_seq.Front();
        dym_val[neuron_model.id_gEInject] += pm.arr_ps[j];  // Add Poisson input
        poisson_time_seq.PopAndFill(pm.arr_pr[j]);
      }
      dbg_printf("Time from %f to %f\n", t_local, t_step_end);
      NextStepSingleNeuronQuiet(dym_val, tmp_neu_state.time_in_refractory[j],
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
    struct TySpikeEvent heading_spike_event(std::numeric_limits<double>::quiet_NaN(), -1);
    struct TyNeuronalDymState<TyNeuronModel> bk_neu_state;
    TySpikeEventVec spike_events;

    dbg_printf("===== NextDt(): t = %f .. %f\n", t, t_end);
    poisson_time_vec.SaveIdxAndClean();
    while (true) {
      spike_events.clear();
      bk_neu_state = neu_state;                    // make a backup
      NextStepNoInteract(neu_state, spike_events, t_end - t);  // Try evolve till the end
      if (spike_events.empty()) {
        t = t_end;
        break;
      } else {
        // Find out the first spike.
        heading_spike_event = *std::min_element(spike_events.begin(), spike_events.end());
        dbg_printf("Neuron [%d] spike at t = %f\n",
                   heading_spike_event.id, heading_spike_event.time);
        // Really evolve the whole system.
        poisson_time_vec.RestoreIdx();
        spike_events.clear();
        neu_state = bk_neu_state;
        NextStepNoInteract(neu_state, spike_events, heading_spike_event.time - t);
        // Record all spike events during the time step above
        // Ideally, there should be only one event: heading_spike_event
        bool b_heading_spike_pushed = false;
        for (size_t i = 0; i < spike_events.size(); i ++) {
          if (spike_events[i].id == heading_spike_event.id) {
            b_heading_spike_pushed = true;
          } else {
            cout << "Unexpected spike (spike before \"first\" spike):  [" << i << "] id = " << spike_events[i].id << " time = " << spike_events[i].time << "\n";
          }
          // Here does not do `ras.emplace_back(spike_events[i])' because
          // effectively the spike interaction are done at the same time.
          ras.emplace_back(heading_spike_event.time, spike_events[i].id);
          vec_n_spike[spike_events[i].id]++;
          SynapticInteraction(neu_state, spike_events[i]);
        }
        // Force the neuron to spike, if not already
        if (!b_heading_spike_pushed) {
          neu_state.dym_vals(heading_spike_event.id, neuron_model.id_V) = neuron_model.Vot_Reset;
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
      if (neu_state.dym_vals(j, neuron_model.id_V) > neuron_model.Vot_Threshold) {
        neu_state.dym_vals(j, neuron_model.id_V) = 0;
        neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();
      }
    }
  }
};

// The event's times for each neuron should be ascending.
void FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path)
{
  std::ifstream fin(path);
  size_t id;
  double time;
  for (auto &i : poisson_time_vec) {
    i.clear();
  }
  while (fin >> id >> time) {
    if (id >= poisson_time_vec.size()) {
      cerr << "FillPoissonEventsFromFile(): number of neuron does not match!" << endl;
      exit(-8);
    } else {
      poisson_time_vec[id].push_back(time);
    }
  }
  for (auto &i : poisson_time_vec) {
    i.push_back(NAN);
  }
  cout << "poisson size[0] = " << poisson_time_vec[0].size() << endl;
  cout << "poisson size[1] = " << poisson_time_vec[1].size() << endl;
}

template<typename TyNeuronModel>
void FillNeuStateFromFile(TyNeuronalDymState<TyNeuronModel> &neu_dym_stat, const char *path)
{
  std::ifstream fin(path);
  double v;
  size_t j = 0;
  neu_dym_stat.dym_vals.setZero();
  while (true) {
    for (int k = 0; fin >> v, k < neu_dym_stat.dym_vals.cols(); k++) {
      neu_dym_stat.dym_vals(j, k) = v;
    }
    if (fin.fail())
      break;
    neu_dym_stat.time_in_refractory[j] = 0;
    j++;
    if (j >= neu_dym_stat.time_in_refractory.size()) {
      cerr << "read " << j << " init data" << endl;
      break;
    }
  }
}

// Read network from text file. The number neurons should be known first.
void FillNetFromPath(TyNeuronalParams &pm, const std::string &name_net)
{
  typedef Eigen::Triplet<double> TyEdgeTriplet;
  std::vector<TyEdgeTriplet> net_coef;
  int n_neu = pm.n_total();

  if (name_net == "-") {
    // Fully connected network
    for (int i = 0; i < n_neu; i++) {
      for (int j = 0; j < n_neu; j++) {
        if (i==j) {
          continue;
        }
        net_coef.push_back(TyEdgeTriplet(i, j, 1.0));
      }
    }
  } else if (name_net == "--") {
    // Simple ring structure
    net_coef.push_back(TyEdgeTriplet(0, n_neu-1, 1.0));
    for (int i = 1; i < n_neu; i++) {
      net_coef.push_back(TyEdgeTriplet(i, i-1, 1.0));
    }
  } else {
    // Read network from text file
    // Negative strength is possible, but that is unusual
    std::ifstream fin_net(name_net);
    double ev;
    for (int i = 0; i < n_neu; i++) {
      for (int j = 0; j < n_neu; j++) {
        fin_net >> ev;
        if (ev != 0 && std::isfinite(ev) && i != j) {
          net_coef.push_back(TyEdgeTriplet(i, j, ev));
        }
      }
    }
  }

  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());
  for (int j = 0; j < pm.n_total(); j++) {
    if (pm.net.coeffRef(j,j)) {
      pm.net.coeffRef(j,j) = 0;
    }
  }
  pm.net.prune(std::numeric_limits<double>::min(), 1);  // remove zeros;
  pm.net.makeCompressed();
}

template<typename TyNeuronModel>
int MainLoop(const po::variables_map &vm)
{
  TyNeuronModel neuron_model;

  // Set neural parameters
  TyNeuronalParams pm(vm["nE"].as<unsigned int>(), vm["nI"].as<unsigned int>());
  int n_neu = vm["nE"].as<unsigned int>() + vm["nI"].as<unsigned int>();

  // set poisson input
  if (! vm.count("ps")) {
    cerr << "Must specify poisson input strength (--ps arg)." << endl;
    return 1;
  }
  for (int i = 0; i < n_neu; i++) {
    pm.arr_pr[i] = vm["pr"].as<double>();
    pm.arr_ps[i] = vm["ps"].as<double>();
  }

//  if (vm.count("pr-mul")) {
//    std::string v = vm["pr-mul"].as<std::string>();
//    cout << v << endl;
//    // TODO: pause it
//  }

  if (vm.count("psi") || vm.count("pri")) {
    cerr << "option --psi and --pri not support yet!" << endl;
  }

  pm.scee = vm["scee"].as<double>();
  pm.scie = vm["scie"].as<double>();
  pm.scei = vm["scei"].as<double>();
  pm.scii = vm["scii"].as<double>();

  if (vm.count("net")) {
    std::string name_net = vm["net"].as<std::string>();
    FillNetFromPath(pm, name_net);
  } else {
    cout << "You must specify the network. (--net)" << endl;
    exit(-1);
  }

  double e_t   = vm["t"].as<double>();
  double e_dt  = vm["dt"].as<double>();
  double e_stv = e_dt;
  if (vm.count("stv")) {
    e_stv = vm["stv"].as<double>();
  }
  if (!(e_t > 0) || !(e_dt > 0) || !(e_stv > 0)) {
    cerr << "Only support positive time!" << endl;
    return 2;
  }
  if (!(e_dt<=e_stv && e_stv<=e_t)) {
    cerr << "Must dt <= stv <= t !" << endl;
    return 2;
  }

  std::ofstream fout_volt;
  bool output_volt = false;
  if (vm.count("volt-path")) {
    output_volt = true;
    fout_volt.open( vm["volt-path"].as<std::string>() );
  }

  std::ofstream fout_G;    // G: conductance
  bool output_G = false;
  if (vm.count("conductance-path")) {
    output_G = true;
    fout_G.open( vm["conductance-path"].as<std::string>() );
  }

  std::ofstream fout_ras;
  bool output_ras = false;
  if (vm.count("ras-path")) {
    output_ras = true;
    fout_ras.open( vm["ras-path"].as<std::string>() );
    fout_ras.precision(17);
  }

  // simulator for the neural network
  NeuronSimulatorExactSpikeSeqEvolve<TyNeuronModel> neu_simu(neuron_model, pm, e_dt);

  if (vm.count("initial-state-path")) {
    FillNeuStateFromFile(neu_simu.neu_state,
                         vm["initial-state-path"].as<std::string>().c_str());
    cout << "initial state loaded!" << endl;
    neu_simu.SaneTestVolt();
  }

  if (vm.count("input-event-path")) {
    FillPoissonEventsFromFile(neu_simu.poisson_time_vec,
                              vm["input-event-path"].as<std::string>().c_str());
    cout << "input event loaded!" << endl;
  }

  if (output_volt) {
    for (int j = 0; j < pm.n_total(); j++) {
      fout_volt.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_V), sizeof(double));
    }
  }
  if (output_G) {
    for (int j = 0; j < pm.n_total(); j++) {
      fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_gE), sizeof(double));
      fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_gI), sizeof(double));
    }
  }

  std::vector<size_t> vec_n_spike(pm.n_total());  // count the number of spikes
  TySpikeEventVec ras;                            // record spike raster
  int n_dt_in_stv = round(e_stv / e_dt);
  int count_n_dt_in_stv = n_dt_in_stv;
  size_t n_step = (size_t)(e_t / e_dt);
  // Main loop
  for (size_t i = 0; i < n_step; i++) {
    neu_simu.NextDt(ras, vec_n_spike);
    if (output_ras) {
      for (size_t j = 0; j < ras.size(); j++) {
        fout_ras << ras[j].id + 1 << '\t' << ras[j].time << '\n';
      }
    }
    ras.clear();
    // output dynamical variable(s) every n_dt_in_stv
    count_n_dt_in_stv--;
    if (count_n_dt_in_stv > 0) {
      continue;
    }
    count_n_dt_in_stv = n_dt_in_stv;
    if (output_volt) {
      for (int j = 0; j < pm.n_total(); j++) {
        fout_volt.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_V), sizeof(double));
      }
    }
    if (output_G) {
      for (int j = 0; j < pm.n_total(); j++) {
        fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_gE), sizeof(double));
        fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, neuron_model.id_gI), sizeof(double));
      }
    }
  }

  if (vm.count("isi-path")) {
    std::ofstream fout_isi( vm["isi-path"].as<std::string>() );
    fout_isi.precision(17);
    for (int j = 0; j < pm.n_total(); j++) {
      fout_isi << e_t / vec_n_spike[j] << '\t';
    }
  }

  return 0;
}

int main(int argc, char *argv[])
{
  // Declare the supported options.
  po::options_description desc("Options");
  // http://stackoverflow.com/questions/3621181/short-options-only-in-boostprogram-options
  desc.add_options()
      ("neuron-model",  po::value<std::string>(),
       "One of LIF-G, LIF-GH")
      ("help,h",
       "produce help message")
      ("t",    po::value<double>()->default_value(1e4),
       "simulation time")
      ("dt",   po::value<double>()->default_value(1.0/2),
       "simulation delta t")
      ("stv",  po::value<double>(),
       "delta t for output")
      ("nE",   po::value<unsigned int>()->default_value(1),
       "number of excitatory neuron")
      ("nI",   po::value<unsigned int>()->default_value(0),
       "number of inhibitory neuron")
      ("net",  po::value<std::string>(),
       "network name")
      ("scee", po::value<double>()->default_value(0.0),
       "cortical strength E->E")
      ("scie", po::value<double>()->default_value(0.0),
       "cortical strength E->I")
      ("scei", po::value<double>()->default_value(0.0),
       "cortical strength I->E")
      ("scii", po::value<double>()->default_value(0.0),
       "cortical strength I->I")
      ("ps",   po::value<double>(),
       "Poisson input strength")
      ("pr",   po::value<double>()->default_value(1.0),
       "Poisson input rate")
      ("psi",  po::value<double>(),
       "Poisson input strength, inhibitory type")
      ("pri",  po::value<double>(),
       "Poisson input rate, inhibitory type")
      ("pr-mul", po::value<double>()->multitoken(),
       "Poisson input rate multiper")
      ("volt-path,o",      po::value<std::string>(),
       "volt output file path")
      ("ras-path",         po::value<std::string>(),
       "ras output file path")
      ("isi-path",         po::value<std::string>(),
       "isi output file path")
      ("conductance-path", po::value<std::string>(),
       "conductance output file path")
      ("initial-state-path", po::value<std::string>(),
       "initial state file path")
      ("input-event-path", po::value<std::string>(),
       "Input event file path")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
  }

  // Set neuron model
  eNeuronModel e_neuron_model = LIF_G;
  if (vm.count("neuron-model")) {
    const std::string &str_nm = vm["neuron-model"].as<std::string>();
    if (str_nm == "LIF-G") {
      e_neuron_model = LIF_G;
    } else if (str_nm == "LIF-GH") {
      e_neuron_model = LIF_GH;
    } else {
      cerr << "Unrecognized neuron model. See help\n";
      return -1;
    }
  } else {
    cout << "Warning: Neuron model not specified, using LIF G model\n";
    e_neuron_model = LIF_G;
  }

  int rt = 0;
  switch (e_neuron_model) {
    case LIF_G:
      rt = MainLoop<Ty_LIF_G>(vm);
      break;
    case LIF_GH:
      rt = MainLoop<Ty_LIF_GH>(vm);
      break;
    default:
      rt = -1;
      cerr << "Unrecognized neuron model enum \n";
  }

  return rt;
}
