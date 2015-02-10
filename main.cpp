#include "common_header.h"
#include "math_helper.h"

// #define NDEBUG
#include <cassert>

using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::SparseMatrix<double> SparseMat;  // column-major by default
typedef std::vector<double> TyArrVals;

std::ofstream fout("a.txt");
//std::ofstream fout("/dev/null");

// list of neuron models
enum NeuronModel
{
  LIF_G,
  LIF_GH,
  HH
};

const struct TyLIFParam
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
  int n_var = 3;  // number of dynamical variables
  int id_V  = 0;  // index of V variable
  int id_gE = 1;  // index of gE variable
  int id_gI = 2;  // index of gI variable
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
  {
    neuron_model = _neuron_model;
    SetNumberOfNeurons(_n_E, _n_I);
  }
};

struct TyNeuronalDymState
{
  // current dynamical states of neurons
  //Eigen::ArrayXXd dym_vals;  // For CNeuronSimulator
  Eigen::Matrix<double, -1, -1, Eigen::RowMajor> dym_vals;  // For CNeuronSimulatorEvolveEach
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

class TyPoissonTimeSeq: public std::vector<double>
{
public:
  iterator it_seq = begin();  // point to current event

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
    it_seq = begin();
  }

  double Front() const
  {
    return *it_seq;
  }

  void PopAndFill(double rate)
  {
    assert(it_seq < end());
    it_seq++;
    if (it_seq == end()) {
      assert(size() > 0);
      if (std::isfinite(back())) {
        Init(rate, back());
        AddEventsUntilTime(rate, back() + 12.0 / rate);
      } else {
        it_seq--;
      }
    }
  }
};

class TyPoissonTimeVec: public std::vector<TyPoissonTimeSeq>
{
public:
  void Init(const TyArrVals &rate_vec, double t0 = 0)
  {
    // each TyPoissonTimeSeq should have at least one event
    resize(rate_vec.size());
    for (size_t j = 0; j < rate_vec.size(); j++) {
      operator[](j).Init(rate_vec[j], t0);
    }
  }
  TyPoissonTimeVec(const TyArrVals &rate_vec, double t0 = 0)
  {
    Init(rate_vec, t0);
  }
  TyPoissonTimeVec() = default;
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

  // evolve assume not reset, not external input at all (isolated)
  void NextDtContinuous(double *dym_val, double &spike_time_local, double dt) const
  {
    // classical Runge Kutta 4 (for voltage)
    // conductance evolve exactly use exp()

    double v_n = dym_val[LIF_param.id_V];
    double k1, k2, k3, k4;
    double exp_E = exp(-0.5 * dt / LIF_param.Time_ExCon);
    double exp_I = exp(-0.5 * dt / LIF_param.Time_InCon);

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
    // k4 = f(t+dt/2, y_n + k3*dt)
    k4 = LIFGetDv(dym_val);

    dym_val[LIF_param.id_V ] = v_n + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);

    if (v_n < LIF_param.Vot_Threshold
        && dym_val[LIF_param.id_V ] >= LIF_param.Vot_Threshold) {
      // could miss some spikes
      spike_time_local = root_search(0, dt,
        v_n, dym_val[LIF_param.id_V ],
        k1, LIFGetDv(dym_val), LIF_param.Vot_Threshold, 1e-12);
    } else {
      spike_time_local = std::numeric_limits<double>::quiet_NaN();
    }
  }

  // evolve condunctance only
  void NextDtConductance(double *dym_val, double dt) const
  {
    dym_val[LIF_param.id_gE] *= exp(-dt / LIF_param.Time_ExCon);
    dym_val[LIF_param.id_gI] *= exp(-dt / LIF_param.Time_InCon);
  }

  void SynapsesInteract(struct TyNeuronalDymState &e_neu_state, const TySpikeEvent &se) const
  {
    if (se.id < pm.n_E) {
      // Excitatory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scee;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scei;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scii;
        }
      }
    }
  }

  // evolve single neuron as if no external input
  void NextQuietDtNeuState(double *dym_val, double &t_in_refractory,
                           double &spike_time_local, double dt_local)
  {
    //! at most one spike allowed during this dt_local
    if (t_in_refractory == 0) {
      NextDtContinuous(dym_val, spike_time_local, dt_local);
      if (!std::isnan(spike_time_local)) {
        t_in_refractory = dt_local - spike_time_local;
        if (t_in_refractory < LIF_param.TIME_REFRACTORY) {
          dym_val[LIF_param.id_V] = LIF_param.VOT_RESET;
        } else {
          // short refractory period, active again
          dt_local = t_in_refractory - LIF_param.TIME_REFRACTORY;
          t_in_refractory = 0;
          // reset dym_val
          NextDtConductance(dym_val, -dt_local);
          dym_val[LIF_param.id_V] = 0;
          double spike_time_local_tmp;
          NextDtContinuous(dym_val, spike_time_local_tmp, dt_local);
          if (!std::isnan(spike_time_local_tmp)) {
            cerr << "NextQuietDtNeuState(): aaaaaaa" << endl;
          }
        }
      }
    } else {
      double t_refractory_remain = LIF_param.TIME_REFRACTORY
                                 - t_in_refractory;
      if (t_refractory_remain < dt_local) {
        NextDtConductance(dym_val, t_refractory_remain);
        // dym_val[LIF_param.id_V] = 0;
        NextDtContinuous(dym_val, spike_time_local, dt_local - t_refractory_remain);
        t_in_refractory = 0;
        if (!std::isnan(spike_time_local)) {
          cerr << "NextQuietDtNeuState(): bbbbbbb" << endl;
        }
      } else {
        NextDtConductance(dym_val, dt_local);
        t_in_refractory += dt_local;
      }
    }
  }

  // evolve as if no synaptic interact
  void NextDt(struct TyNeuronalDymState &tmp_neu_state,
              TySpikeEventQueue &spike_events, double dt)
  {
    double t_step_end = t + dt;
    // evolve each neuron as if no reset
    for (int j = 0; j < pm.n_total(); j++) {
      //! tmp_neu_state.dym_vals must be Row major !
      //TyPoissonTimeQueue &poisson_time_seq = poisson_time_vec[j];
      TyPoissonTimeSeq poisson_time_seq = poisson_time_vec[j];  // an copy
      double *dym_val = tmp_neu_state.dym_vals.data() + j * LIF_param.n_var;
      double t_local = t;
      double spike_time_local = std::numeric_limits<double>::quiet_NaN();
      while (poisson_time_seq.Front() < t_step_end) {
        double dt_local = poisson_time_seq.Front() - t_local;
        NextQuietDtNeuState(dym_val, tmp_neu_state.time_in_refractory[j],
                            spike_time_local, dt_local);
        if (!std::isnan(spike_time_local)) {
          spike_events.emplace(t_local + spike_time_local, j);
        }
        dym_val[LIF_param.id_gE] += pm.arr_ps[j];
        poisson_time_seq.PopAndFill(pm.arr_pr[j]);
        t_local += dt_local;
      }
      NextQuietDtNeuState(dym_val, tmp_neu_state.time_in_refractory[j],
                          spike_time_local, t_step_end - t_local);
      if (!std::isnan(spike_time_local)) {
        spike_events.emplace(t_local + spike_time_local, j);
      }
    }
  }

  // evolve the whole system one dt
  void NextStep()
  {
    double t_sub = t;
    double t_end = t + dt;
    TySpikeEvent latest_spike_event =
      {std::numeric_limits<double>::quiet_NaN(), -1};
    struct TyNeuronalDymState tmp_neu_state;

    do {
      cout << "t = " << t_sub << endl;
      TySpikeEventQueue spike_events;
      tmp_neu_state = neu_state;
      NextDt(tmp_neu_state, spike_events, t_end - t_sub);
      if (spike_events.size() == 0) {
        neu_state = tmp_neu_state;
        break;
      } else {
        cout << "spike [" << spike_events.top().id << "] at "
             << spike_events.top().time << endl;
        latest_spike_event = spike_events.top();
        // Really evolve the whole system, spike_events will be discard.
        NextDt(neu_state, spike_events, latest_spike_event.time - t_sub);
        if (spike_events.top().id != latest_spike_event.id) {
          cerr << "NextStep(): Should not happen." << endl;
        }
        t_sub = latest_spike_event.time;
        // force the neuron to spike
        if (neu_state.time_in_refractory[latest_spike_event.id] == 0) {
          neu_state.dym_vals(latest_spike_event.id, LIF_param.id_V) = 0;
          neu_state.time_in_refractory[latest_spike_event.id] =
            std::numeric_limits<double>::min();
        }
        SynapsesInteract(neu_state, latest_spike_event);
      }
    } while (true);
    t = t_end;
  }
};

typedef CNeuronSimulatorEvolveEach CNeuronSimulator;

int main()
{
  // fill in parameters
  int n_neu = 2;
  TyNeuronalParams pm(LIF_G, n_neu, 0);
  for (int i = 0; i < n_neu; i++) {
    pm.arr_pr[i] = 1.0;
    pm.arr_ps[i] = 0.012;
  }
  std::vector< Eigen::Triplet<double> > net_coef;
  //net_coef.push_back({1, 0, 1.0});
  for (int i = 0; i < n_neu; i++) {
    for (int j = 0; j < n_neu; j++) {
      if (i==j) {
        continue;
      }
      net_coef.push_back({i, j, 1.0});
    }
  }
  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());
  for (int j = 0; j < pm.n_total(); j++) {
    if (pm.net.coeffRef(j,j)) {
      pm.net.coeffRef(j,j) = 0;
    }
  }
  pm.net.makeCompressed();

  // Create simulator
  double e_dt = 1.0 / 2;  // ms
  CNeuronSimulator neu_simu(pm, e_dt);

  cout << "t = " << neu_simu.t << endl;
  cout << neu_simu.neu_state.dym_vals << endl;

//  for (int i = 0; i < (int)(1e3 / e_dt); i++) {
//    neu_simu.NextStep();
//  }

  for (int i = 0; i < 50; i++) {
    neu_simu.NextStep();
    cout << "t = " << neu_simu.t << endl;
    cout << neu_simu.neu_state.dym_vals << endl;
  }

  // Octave
  // load a.txt
  // plot([a(1:2:end, 1), a(2:2:end, 1)], '-o')

  return 0;
}
