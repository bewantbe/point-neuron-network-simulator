#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using std::cout;
using std::cerr;
using std::endl;

// list of neuron models
enum NeuronModel
{
  LIF_G,
  LIF_GH,
  HH
};

typedef Eigen::SparseMatrix<double> SparseMat;
typedef std::vector<double> TyArrVals;

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
  Eigen::ArrayXXd dym_vals;  // current dynamical states of neurons
  TyArrVals time_in_refractory;

  TyNeuronalDymState(const TyNeuronalParams &pm)
  {
    switch (pm.neuron_model) {
      case LIF_G:
        dym_vals.resize(pm.n_total(), 3);  // V, G_E, G_I
        break;
      case LIF_GH:
        break;
      case HH:
        break;
      default:
        cerr << "Neuron model not support!" << endl;
        exit(-1);
    }
    time_in_refractory.resize(pm.n_total());
  }

//private:
  TyNeuronalDymState()  // prevent empty initialization
  { }
};

//std::random_device rd;
std::mt19937 rand_eng(1);  // gen(rd())

//class TyPoissonQueue
//{
//public:
//  std::queue<double> poisson_time;
//  void AddEventsUntilTime(double rate, double t_until)
//  {
//    std::exponential_distribution<> exp_dis(rate);
//    while (poisson_time.back() <= t_until) {
//      poisson_time.push( poisson_time.back() + exp_dis(rand_eng) );
//    }
//  }
//
//  void Init(double rate, double t0)
//  {
//    std::exponential_distribution<> exp_dis(rate);
//    poisson_time.push(t0 + exp_dis(rand_eng));
//  }
//};

//class TyPoissonQueueVec
//{
//public:
//  typedef std::pair<double, int> TySpikeEvent;  // event time, index (of neuron)
//  struct std::vector<TyPoissonQueue> poission_queue_vec;
//  struct std::vector<TyPoissonQueue>::iterator poission_queue_it;
//  std::vector<TySpikeEvent> poisson_event_vec;
//
//  void Init(const TyArrVals &rate_vec, double t0 = 0)
//  {
//    // each poission_queue_vec[] should have at least one event
//    poission_queue_vec.resize(rate_vec.size());
//    for (size_t j = 0; j < rate_vec.size(); j++) {
//      poission_queue_vec[j].Init(rate_vec[j], t0);
//    }
//    poission_queue_it = poission_queue_vec.end();  // point to next event
//  }
//
//  TySpikeEvent & PopFirst(const TyArrVals &rate_vec)
//  {
//    if (poission_queue_it == poisson_event_vec.end()) {
//      // no enough event, generate more
//      for (int j = 0; j < rate_vec.size(); j++) {
//        poission_queue_vec[j].AddEventsUntilTime(rate_vec[j], 100);
//      }
//    } else {
//      return *poisson_queue_it++;
//    }
//  }
//};

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
  {}

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

class CPoissonEvents
{
public:
  std::priority_queue<
    TySpikeEvent,
    std::vector<TySpikeEvent>,
    std::greater<TySpikeEvent> > spike_events;

  void Clear()
  {
    while (!spike_events.empty()) {
      spike_events.pop();
    }
  }

  void Init(const TyArrVals &rate_vec, double t0)
  {
    Clear();
    for (size_t j = 0; j < rate_vec.size(); j++) {
      if (rate_vec[j]==0) {
        spike_events.emplace(INFINITY, j);
      } else if (rate_vec[j]<0) {
        spike_events.emplace(INFINITY, j);  // anyway ...
        std::cerr << "Negative Poisson rate!" << endl;
      } else {
        // "1.0 - " to avoid log(zero)
        spike_events.emplace(j,
          t0 - log(1.0 - g_rand()) / rate_vec[j]);
      }
    }
  }

  const TySpikeEvent & HeadEvent() const
  {
    return spike_events.top();
  }

  void PopEvent(const TyArrVals &rate_vec)
  {
    const TySpikeEvent event = spike_events.top();
    spike_events.pop();
    spike_events.emplace(event.id,
      event.time - log(1.0 - g_rand()) / rate_vec[event.id]);
  }

  CPoissonEvents(const TyArrVals &rate_vec, double t0  = 0)
  {
    Init(rate_vec, t0);
  }

  CPoissonEvents() {}
};

class CNeuronSimulator
{
public:
  struct TyNeuronalParams pm;
  struct TyNeuronalDymState neu_state;
  CPoissonEvents poisson_events;
  double t, dt;

  CNeuronSimulator(const TyNeuronalParams &_pm, double _dt)
  :pm(_pm), neu_state(_pm)
  {
    t = 0;
    dt = _dt;
    poisson_events.Init(pm.arr_pr, t);
  }

  struct TyNeuronalDymState tmp_neu_state;

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
  } LIF_param;

  // Get instantaneous dv/dt for LIF model
  void LIFGetDvVec(Eigen::ArrayXd &dv_vec, const struct TyNeuronalDymState &nds)
  {
    const Eigen::ArrayXd &v  = nds.dym_vals.col(0);
    const Eigen::ArrayXd &gE = nds.dym_vals.col(1);
    const Eigen::ArrayXd &gI = nds.dym_vals.col(2);

    dv_vec = - LIF_param.Con_Leakage * (v - LIF_param.Vot_Leakage)
             - gE * (v - LIF_param.Vot_Excitatory)
             - gI * (v - LIF_param.Vot_Inhibitory);
  }

  // evolve assume not reset, not external input at all
  void NextDtContinuous(struct TyNeuronalDymState &nds, double dt)
  {
    // classical Runge Kutta 4
    // for voltage only, conductance evolve exactly use exp()

    struct TyNeuronalDymState tmp_nds = nds;
    Eigen::ArrayXd k1(pm.n_total());
    Eigen::ArrayXd k2(pm.n_total());
    Eigen::ArrayXd k3(pm.n_total());
    Eigen::ArrayXd k4(pm.n_total());

    // k1 = f(t_n, y_n)
    LIFGetDvVec(k1, tmp_nds);

    // y_n + 0.5*k1*dt
    tmp_nds.dym_vals.col(0) = nds.dym_vals.col(0) + 0.5 * dt * k1;
    tmp_nds.dym_vals.col(1) = nds.dym_vals.col(1) * exp(-0.5 * dt / LIF_param.Time_ExCon);
    tmp_nds.dym_vals.col(2) = nds.dym_vals.col(2) * exp(-0.5 * dt / LIF_param.Time_InCon);
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    LIFGetDvVec(k2, tmp_nds);

    // y_n + 0.5*k2*dt
    tmp_nds.dym_vals.col(0) = nds.dym_vals.col(0) + 0.5 * dt * k2;
    tmp_nds.dym_vals.col(1) = nds.dym_vals.col(1) * exp(-0.5 * dt / LIF_param.Time_ExCon);
    tmp_nds.dym_vals.col(2) = nds.dym_vals.col(2) * exp(-0.5 * dt / LIF_param.Time_InCon);
    // k3 = f(t+dt/2, y_n + 0.5*k2*dt)
    LIFGetDvVec(k3, tmp_nds);

    // y_n + 0.5*k3*dt
    tmp_nds.dym_vals.col(0) = nds.dym_vals.col(0) + dt * k3;
    tmp_nds.dym_vals.col(1) = nds.dym_vals.col(1) * exp(- dt / LIF_param.Time_ExCon);
    tmp_nds.dym_vals.col(2) = nds.dym_vals.col(2) * exp(- dt / LIF_param.Time_InCon);
    // k4 = f(t+dt, y_n + k3*dt)
    LIFGetDvVec(k4, tmp_nds);

    nds.dym_vals.col(0) += dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);
    nds.dym_vals.col(1) = tmp_nds.dym_vals.col(1);
    nds.dym_vals.col(2) = tmp_nds.dym_vals.col(2);

    t += dt;
  }

  void NextDt(double dt)
  {
    tmp_neu_state = neu_state;
    // evolve as if not reset not resting
    NextDtContinuous(tmp_neu_state, dt);
    // keep fired neuron in resting state;
    for (int j = 0; j < pm.n_total(); j++) {
      // TODO: not exactly correct
      if (tmp_neu_state.time_in_refractory[j]) {
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    // do reset
    for (int j = 0; j < pm.n_total(); j++) {
      // TODO: not exactly correct
      if (tmp_neu_state.dym_vals(j, 0) > LIF_param.Vot_Threshold) {
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    neu_state = tmp_neu_state;
  }

  // evolve one dt
  void NextStep()
  {
    double t_step_end = t + dt;
    TySpikeEvent se;
    while ((se = poisson_events.HeadEvent()).time < t_step_end) {
      NextDt(se.time - t);        // step one sub dt
      poisson_events.PopEvent(pm.arr_pr);
    }
  }
};

int main()
{
  // fill in parameters
  TyNeuronalParams pm(LIF_G, 2, 0);
  pm.arr_pr[0] = 1.0;
  pm.arr_pr[1] = 1.0;
  pm.arr_ps[0] = 0.012;
  pm.arr_ps[1] = 0.012;
  std::vector< Eigen::Triplet<double> > net_coef;
  net_coef.push_back({1, 0, 1.0});
  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());

  // Create simulator
  double e_dt = 1.0 / 2;  // ms
  CNeuronSimulator neu_simu(pm, e_dt);

  for (int i = 0; i < 1; i++) {
    neu_simu.NextStep();
    cout << "t = " << neu_simu.t << endl;
    cout << neu_simu.neu_state.dym_vals << endl;
  }

  return 0;
}
