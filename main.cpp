// g++ -g -O0 -std=c++11 -Wall main.cpp math_helper.cpp -lboost_program_options -o bin/vec_IFsimu
//Could faster:
// g++ -g -O2 -falign-functions=16 -falign-loops=16
#include "common_header.h"
#include "math_helper.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <cassert>

const bool g_b_debug = false;
#define NDEBUG
#ifdef NDEBUG
#define dbg_printf(...) ((void)0);
#else
#define dbg_printf printf
#endif

#include "fmath.hpp"

using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMat;
typedef std::vector<double> TyArrVals;

// list of neuron models
enum NeuronModel
{
  LIF_G,
  LIF_GH,
  HH
};

struct TyLIFParam
{
  double Vot_Threshold  = 1.0;
  double VOT_RESET      = 0.0;
  double Con_Leakage = 0.05;
  double Vot_Leakage = 0.0;
  double Vot_Excitatory = 14.0/3.0;
  double Vot_Inhibitory = -2.0/3.0;
  double Time_ExCon    = 2.0;
  double Time_ExConR   = 0.5;
  double Time_InCon    = 5.0;
  double Time_InConR   = 0.8;
  double TIME_REFRACTORY = 2.0;  // ms
  static const int n_var = 3;  // number of dynamical variables
  static const int id_V  = 0;  // index of V variable
  static const int id_gE = 1;  // index of gE variable
  static const int id_gI = 2;  // index of gI variable
  static const int id_gEInject = 1;  // index of gE injection variable
  static const int id_gIInject = 2;  // index of gI injection variable
} LIF_param;

// parameters about neuronal system
struct TyNeuronalParams
{
  enum NeuronModel neuron_model;
  int n_E, n_I;
  SparseMat net;
  double scee, scie, scei, scii;
  TyArrVals arr_pr;
  TyArrVals arr_ps;
  TyArrVals arr_psi;

  int n_total() const
  { return n_E + n_I; }

  // should use this this function to initialize neurons
  void SetNumberOfNeurons(int _n_E, int _n_I)
  {
    n_E = _n_E;
    n_I = _n_I;
    arr_pr.resize(n_total());
    arr_ps.resize(n_total());
    arr_psi.resize(n_total());
    net.resize(n_total(), n_total());
  }

  TyNeuronalParams(enum NeuronModel _neuron_model, int _n_E, int _n_I)
  :scee(0), scie(0), scei(0), scii(0)
  {
    neuron_model = _neuron_model;
    SetNumberOfNeurons(_n_E, _n_I);
  }
};

struct TyNeuronalDymState
{
  // current dynamical states of neurons
  // RowMajor for CNeuronSimulatorEvolveEach, so state variables
  // of each neuron are continuous on memory
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dym_vals;
  TyArrVals time_in_refractory;

  void Zeros()
  {
    dym_vals.setZero();
    for (auto &i : time_in_refractory) i = 0;
  }

  TyNeuronalDymState(const TyNeuronalParams &pm)
  {
    switch (pm.neuron_model) {
      case LIF_G:
        dym_vals.resize(pm.n_total(), LIF_param.n_var);  // V, G_E, G_I
        break;
      case LIF_GH:
//        break;
      case HH:
//        break;
      default:
        cerr << "Neuron model not support!" << endl;
        exit(-1);
    }
    time_in_refractory.resize(pm.n_total());
    Zeros();
  }

//private:
  TyNeuronalDymState()  // prevent empty initialization
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

class TyPoissonTimeSeq: public TyInternalVec
{
public:
  typedef size_t TyIdx;
  TyIdx id_seq = 0;  // point to current event
  TyPoissonTimeSeq() = default;
  TyPoissonTimeSeq(const TyPoissonTimeSeq &pts)
    : TyInternalVec(pts),
    id_seq(pts.id_seq)
  {}

  // Requirements: t_until <= inf, rate >= 0.
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

class TyPoissonTimeVec: public std::vector<TyPoissonTimeSeq>
{
  std::vector<TyPoissonTimeSeq::TyIdx> id_seq_vec;
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
      if (it->size() - it->id_seq < it->id_seq / 7) {
        it->Shrink();
        dbg_printf("  ( %d-th TyPoissonTimeSeq size shrinked )\n", j);
      }
      id_seq_vec[j] = it->id_seq;
    }
  }
};

class CNeuronSimulatorEvolveEach
{
public:
  struct TyNeuronalParams pm;
  struct TyNeuronalDymState neu_state;
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  CNeuronSimulatorEvolveEach(const TyNeuronalParams &_pm, double _dt)
  :pm(_pm), neu_state(_pm)
  {
    dt = _dt;
    t = 0;
    poisson_time_vec.Init(pm.arr_pr, t);
  }

  // Get instantaneous dv/dt for LIF model
  double LIFGetDv(const double *dym_val) const
  {
    return - LIF_param.Con_Leakage * (dym_val[LIF_param.id_V] - LIF_param.Vot_Leakage)
           - dym_val[LIF_param.id_gE] * (dym_val[LIF_param.id_V] - LIF_param.Vot_Excitatory)
           - dym_val[LIF_param.id_gI] * (dym_val[LIF_param.id_V] - LIF_param.Vot_Inhibitory);
  }

  // Evolve the ODE assume not reset, not external input at all.
  // Return derivative k1 at t_n, for later interpolation.
  __attribute__ ((noinline)) double LIFDymRK4Step(double *dym_val, double dt) const
  {
    // classical Runge Kutta 4 (for voltage)
    // conductance evolve exactly use exp()

    double v_n = dym_val[LIF_param.id_V];
    double k1, k2, k3, k4;
    double exp_E = fmath::expd(-0.5 * dt / LIF_param.Time_ExCon);
    double exp_I = fmath::expd(-0.5 * dt / LIF_param.Time_InCon);

    // k1 = f(t_n, y_n)
    k1 = LIFGetDv(dym_val);

    // y_n + 0.5*k1*dt
    dym_val[LIF_param.id_V ] = v_n + 0.5 * dt * k1;
    dym_val[LIF_param.id_gE] *= exp_E;
    dym_val[LIF_param.id_gI] *= exp_I;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k2 = LIFGetDv(dym_val);

    // y_n + 0.5*k2*dt
    dym_val[LIF_param.id_V ] = v_n + 0.5 * dt * k2;
    // k3 = f(t+dt/2, y_n + 0.5*k2*dt)
    k3 = LIFGetDv(dym_val);

    // y_n + k3*dt
    dym_val[LIF_param.id_V ] = v_n + dt * k3;
    dym_val[LIF_param.id_gE] *= exp_E;
    dym_val[LIF_param.id_gI] *= exp_I;
    // k4 = f(t+dt, y_n + k3*dt)
    k4 = LIFGetDv(dym_val);

    dym_val[LIF_param.id_V ] = v_n + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);

    return k1;
  }

  // Evolve the ODE assume not reset, not external input at all.
  // With spike time calculation.
  // `spike_time_local' should be guaranteed to be with in [0, dt] or NAN.
  // The separation of LIFDymRK4Step() make program 4% slower.
  __attribute__ ((noinline)) void NextDtContinuous(double *dym_val, double &spike_time_local, double dt) const
  {
    double dym_tn[LIF_param.n_var];
    memcpy(dym_tn, dym_val, LIF_param.n_var*sizeof(double));

    double v_n = dym_val[LIF_param.id_V];
    double k1  = LIFDymRK4Step(dym_val, dt);

    if (v_n <= LIF_param.Vot_Threshold
        && dym_val[LIF_param.id_V ] > LIF_param.Vot_Threshold) {
      spike_time_local = cubic_hermit_real_root(dt,
        v_n, dym_val[LIF_param.id_V ],
        k1, LIFGetDv(dym_val), LIF_param.Vot_Threshold);
//      // Further root refinement
//      double dym_tmp[LIF_param.n_var];
//      memcpy(dym_tmp, dym_tn, LIF_param.n_var*sizeof(double));
//      LIFDymRK4Step(dym_tmp, spike_time_local);
//      double dr = (dym_tmp[LIF_param.id_V] - LIF_param.Vot_Threshold) / LIFGetDv(dym_tmp);
//      static double dr_max = 0;
//      if (fabs(dr) > 0.0005) {
//        dr_max = dr;
//        cout.precision(16);
//        cout << "dr_max = " << dr_max << ", t = " << t << "\n";
//        cout << "  dt  = " << dt << "\n";
//        cout << "  @t_n  : v, gE, gI = " << dym_tn[0] << " " << dym_tn[1]<< " " << dym_tn[2] << "\n";
//        cout << "  @t_n+1: v, gE, gI = " << dym_val[0] << " " << dym_val[1]<< " " << dym_val[2] << "\n";
//        cout << "  k_n   = " << LIFGetDv(dym_tn) << "\n";
//        cout << "  k_n+1 = " << LIFGetDv(dym_val) << "\n";
//        cout << "  spike_time_local = " << spike_time_local << "\n";
//      }
//      spike_time_local -= dr;
//      throw "whhhhhhh";
    } else {
      if (v_n > 0.996 && k1>0) { // the v_n > 0.996 is for dt=0.5 ms, LIF,G model
        // Try capture some missing spikes that the intermediate value passes
        // threshold, but both ends lower than threshold.
        // Get quadratic curve from value of both ends and derivative from left end
        // Return the maximum point as `t_max_guess'
        double c = v_n;
        double b = k1;
        double a = (dym_val[LIF_param.id_V ] - c - b*dt)/(dt*dt);
        double t_max_guess = -b/(2*a);
        // In LIF with jump conductance, it can guarantee that a<0 (concave),
        // hence t_max_guess > 0. But with G H conductance, we still need to
        // check 0 < t_max_guess
        if (0 < t_max_guess && t_max_guess<dt
            && (b*b)/(-4*a)+c > LIF_param.Vot_Threshold) {
          //dbg_printf("Rare event: mid-dt spike captured, guess time: %f\n", t_max_guess);
          printf("NextDtContinuous(): possible mid-dt spike detected:\n");
          printf("  Guessed max time: %f, t = %f, dt = %f\n", t_max_guess, t, dt);
          // root should in [0, t_max_guess]
          spike_time_local = cubic_hermit_real_root(dt,
            v_n, dym_val[LIF_param.id_V ],
            k1, LIFGetDv(dym_val), LIF_param.Vot_Threshold);
        } else {
          spike_time_local = std::numeric_limits<double>::quiet_NaN();
        }
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // evolve condunctance only
  void NextDtConductance(double *dym_val, double dt) const
  {
    dym_val[LIF_param.id_gE] *= fmath::expd(-dt / LIF_param.Time_ExCon);
    dym_val[LIF_param.id_gI] *= fmath::expd(-dt / LIF_param.Time_InCon);
  }

  __attribute__ ((noinline)) void SynapticInteraction(struct TyNeuronalDymState &e_neu_state, const TySpikeEvent &se) const
  {
    if (se.id < pm.n_E) {
      // Excitatory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gEInject) += it.value() * pm.scee;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gEInject) += it.value() * pm.scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gIInject) += it.value() * pm.scei;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gIInject) += it.value() * pm.scii;
        }
      }
    }
  }

  // Evolve single neuron as if no external input.
  // Return first spike time in `spike_time_local', if any.
  __attribute__ ((noinline)) void NextQuietDtNeuState(double *dym_val, double &t_in_refractory,
                           double &spike_time_local, double dt_local)
  {
    //! at most one spike allowed during this dt_local
    if (t_in_refractory == 0) {
      dbg_printf("NextQuietDtNeuState(): dt_local = %f:\n", dt_local);
      dbg_printf("NextQuietDtNeuState(): begin state=%f,%f,%f\n",
                 dym_val[0], dym_val[1], dym_val[2]);
      NextDtContinuous(dym_val, spike_time_local, dt_local);
      dbg_printf("NextQuietDtNeuState(): end   state=%f,%f,%f\n",
                 dym_val[0], dym_val[1], dym_val[2]);
      if (!std::isnan(spike_time_local)) {
        // Add `numeric_limits<double>::min()' to make sure t_in_refractory > 0.
        t_in_refractory = dt_local - spike_time_local
                          + std::numeric_limits<double>::min();
        dym_val[LIF_param.id_V] = LIF_param.VOT_RESET;
        dbg_printf("Reach threshold detected\n");
        if (t_in_refractory >= LIF_param.TIME_REFRACTORY) {
          // Short refractory period (< dt_local), neuron will active again.
          dt_local = t_in_refractory - LIF_param.TIME_REFRACTORY;
          t_in_refractory = 0;
          // Back to the activation time.
          NextDtConductance(dym_val, -dt_local);
          double spike_time_local_tmp;
          NextDtContinuous(dym_val, spike_time_local_tmp, dt_local);
          if (!std::isnan(spike_time_local_tmp)) {
            cerr << "NextQuietDtNeuState(): Multiple spikes in one dt. Interaction dropped." << endl;
            cerr << "  dt_local = " << dt_local << '\n';
            cerr << "  spike_time_local = " << spike_time_local << '\n';
            cerr << "  t_in_refractory = " << t_in_refractory << '\n';
            cerr << "  dym_val[LIF_param.id_V] = " << dym_val[LIF_param.id_V] << '\n';
            cerr << "  dym_val[LIF_param.id_gE] = " << dym_val[LIF_param.id_gE] << '\n';
            throw "Multiple spikes in one dt.";
          }
        }
      }
    } else {
      // Neuron in refractory period.
      double dt_refractory_remain = LIF_param.TIME_REFRACTORY
                                 - t_in_refractory;
      if (dt_refractory_remain < dt_local) {
        // neuron will awake after dt_refractory_remain which is in this dt_local
        NextDtConductance(dym_val, dt_refractory_remain);
        assert( dym_val[LIF_param.id_V] == LIF_param.VOT_RESET );
        t_in_refractory = 0;
        NextQuietDtNeuState(dym_val, t_in_refractory,
                            spike_time_local, dt_local - dt_refractory_remain);
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
        NextDtConductance(dym_val, dt_local);
        t_in_refractory += dt_local;
      }
    }
  }

  // Evolve every neurons as if no synaptic interaction
  __attribute__ ((noinline)) void NextDt(struct TyNeuronalDymState &tmp_neu_state,
              TySpikeEventVec &spike_events, double dt)
  {
    double t_step_end = t + dt;
    dbg_printf("----- Trying t = %f .. %f -------\n", t, t_step_end);
    for (int j = 0; j < pm.n_total(); j++) {
      //! tmp_neu_state.dym_vals must be Row major !
      TyPoissonTimeSeq &poisson_time_seq = poisson_time_vec[j];
      double *dym_val = tmp_neu_state.dym_vals.data() + j * LIF_param.n_var;
      double t_local = t;
      double spike_time_local = std::numeric_limits<double>::quiet_NaN();
      while (poisson_time_seq.Front() < t_step_end) {
        dbg_printf("Neuron %d receive %lu-th (%lu) Poisson input at %f\n",
                   j, poisson_time_seq.id_seq, poisson_time_seq.size(),
                   poisson_time_seq.Front());
        dbg_printf("Time from %f to %f\n", t_local, poisson_time_seq.Front());
        NextQuietDtNeuState(dym_val, tmp_neu_state.time_in_refractory[j],
                            spike_time_local, poisson_time_seq.Front() - t_local);
        if (!std::isnan(spike_time_local)) {
          spike_events.emplace_back(t_local + spike_time_local, j);
        }
        t_local = poisson_time_seq.Front();
        dym_val[LIF_param.id_gEInject] += pm.arr_ps[j];  // Add Poisson input
        poisson_time_seq.PopAndFill(pm.arr_pr[j]);
      }
      dbg_printf("Time from %f to %f\n", t_local, t_step_end);
      NextQuietDtNeuState(dym_val, tmp_neu_state.time_in_refractory[j],
                          spike_time_local, t_step_end - t_local);
      if (!std::isnan(spike_time_local)) {
        spike_events.emplace_back(t_local + spike_time_local, j);
      }
      dbg_printf("Neuron %d @t=%f end state %f, %f, %f\n", j, t_step_end, dym_val[0], dym_val[1], dym_val[2]);
    }
  }

  // Evolve the whole system one dt, with interaction
  __attribute__ ((noinline)) void NextStep(TySpikeEventVec &ras, std::vector< size_t > &vec_n_spike)
  {
    double t_end = t + dt;
    TySpikeEvent heading_spike_event =
      {std::numeric_limits<double>::quiet_NaN(), -1};
    struct TyNeuronalDymState bk_neu_state;
    TySpikeEventVec spike_events;

    dbg_printf("===== NextStep(): t = %f .. %f\n", t, t_end);
    poisson_time_vec.SaveIdxAndClean();
    while (true) {
      spike_events.clear();
      bk_neu_state = neu_state;                    // make a backup
      NextDt(neu_state, spike_events, t_end - t);  // Try evolve till the end
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
        NextDt(neu_state, spike_events, heading_spike_event.time - t);
        bool b_heading_spike_pushed = false;
        for (size_t i = 0; i < spike_events.size(); i ++) {
          if (spike_events[i].id == heading_spike_event.id) {
            b_heading_spike_pushed = true;
          } else {
            cout << "Non-heading spike:  [" << i << "] id = " << spike_events[i].id << " time = " << spike_events[i].time << "\n";
          }
          // Here does not do `ras.emplace_back(spike_events[i])' because
          // effectively the spike interaction are done at the same time.
          ras.emplace_back(heading_spike_event.time, spike_events[i].id);
          vec_n_spike[spike_events[i].id]++;
          SynapticInteraction(neu_state, spike_events[i]);
        }
        // force the neuron to spike, if not already
        if (!b_heading_spike_pushed) {
          neu_state.dym_vals(heading_spike_event.id, LIF_param.id_V) = LIF_param.VOT_RESET;
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
      if (neu_state.dym_vals(j, LIF_param.id_V) > LIF_param.Vot_Threshold) {
        neu_state.dym_vals(j, LIF_param.id_V) = 0;
        neu_state.time_in_refractory[j] = std::numeric_limits<double>::min();
      }
    }
  }
};

typedef CNeuronSimulatorEvolveEach CNeuronSimulator;

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

void FillNeuStateFromFile(TyNeuronalDymState &neu_dym_stat, const char *path)
{
  std::ifstream fin(path);
  double v, gE, gI;
  size_t j = 0;
  neu_dym_stat.dym_vals.setZero();
  while (fin >> v >> gE >> gI) {
    neu_dym_stat.dym_vals(j, 0) = v;
    neu_dym_stat.dym_vals(j, 1) = gE;
    neu_dym_stat.dym_vals(j, 2) = gI;
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

int main(int argc, char *argv[])
{
  // Declare the supported options.
  po::options_description desc("Options");
  // http://stackoverflow.com/questions/3621181/short-options-only-in-boostprogram-options
  desc.add_options()
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

  // Set neural parameters
  TyNeuronalParams pm(LIF_G, vm["nE"].as<unsigned int>(), vm["nI"].as<unsigned int>());
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
  CNeuronSimulator neu_simu(pm, e_dt);

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
      fout_volt.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_V), sizeof(double));
    }
  }
  if (output_G) {
    for (int j = 0; j < pm.n_total(); j++) {
      fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_gE), sizeof(double));
      fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_gI), sizeof(double));
    }
  }

  std::vector<size_t> vec_n_spike(pm.n_total());  // count the number of spikes
  TySpikeEventVec ras;                            // record spike raster
  int n_dt_in_stv = round(e_stv / e_dt);
  int count_n_dt_in_stv = n_dt_in_stv;
  size_t n_step = (size_t)(e_t / e_dt);
  // Main loop
  for (size_t i = 0; i < n_step; i++) {
    neu_simu.NextStep(ras, vec_n_spike);
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
        fout_volt.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_V), sizeof(double));
      }
    }
    if (output_G) {
      for (int j = 0; j < pm.n_total(); j++) {
        fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_gE), sizeof(double));
        fout_G.write((char*)&neu_simu.neu_state.dym_vals(j, LIF_param.id_gI), sizeof(double));
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
