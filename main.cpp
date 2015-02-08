#include "common_header.h"

using std::cout;
using std::cerr;
using std::endl;

std::ofstream fout("a.txt");
//std::ofstream fout("/dev/null");

// list of neuron models
enum NeuronModel
{
  LIF_G,
  LIF_GH,
  HH
};

typedef Eigen::SparseMatrix<double> SparseMat;  // column-major by default
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

struct TyNeuronalDymState
{
  Eigen::ArrayXXd dym_vals;  // current dynamical states of neurons
  TyArrVals time_in_refractory;

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

typedef std::priority_queue<
    TySpikeEvent,
    std::vector<TySpikeEvent>,
    std::greater<TySpikeEvent> > TySpikeEventQueue;

class CPoissonEvents
{
public:
  TySpikeEventQueue poisson_spike_events;

  void Clear()
  {
    while (!poisson_spike_events.empty()) {
      poisson_spike_events.pop();
    }
  }

  void Init(const TyArrVals &rate_vec, double t0)
  {
    Clear();
    for (size_t j = 0; j < rate_vec.size(); j++) {
      if (rate_vec[j]==0) {
        poisson_spike_events.emplace(INFINITY, j);
      } else if (rate_vec[j]<0) {
        poisson_spike_events.emplace(INFINITY, j);  // anyway ...
        std::cerr << "Negative Poisson rate!" << endl;
      } else {
        // "1.0 - " to avoid log(zero)
        poisson_spike_events.emplace(
          t0 - log(1.0 - g_rand()) / rate_vec[j],
          j);
      }
    }
  }

  const TySpikeEvent & HeadEvent() const
  {
    return poisson_spike_events.top();
  }

  void PopEvent(const TyArrVals &rate_vec)
  {
    const TySpikeEvent event = poisson_spike_events.top();
    poisson_spike_events.pop();
    poisson_spike_events.emplace(
      event.time - log(1.0 - g_rand()) / rate_vec[event.id],
      event.id);
  }

  CPoissonEvents(const TyArrVals &rate_vec, double t0  = 0)
  {
    Init(rate_vec, t0);
  }

  CPoissonEvents() {}
};

//class TyTimeQueue
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

class TyTimeQueue: private std::queue<double>
{
public:
  void AddEventsUntilTime(double rate, double t_until)
  {
    std::exponential_distribution<> exp_dis(rate);
    while (back() <= t_until) {
      push( back() + exp_dis(rand_eng) );
    }
  }

  void Init(double rate, double t0)
  {
    std::exponential_distribution<> exp_dis(rate);
    push(t0 + exp_dis(rand_eng));
  }

  double Front() const
  {
    return front();
  }

  void PopWithRateTime(double rate)
  {
    if (size() <= 1) {
      if (std::isfinite(back())) {
        AddEventsUntilTime(rate, back() + 10.0 / rate);
        pop();
      }
    } else {
      pop();
    }
  }
};

class TyPoissonTimeVec: public std::vector<TyTimeQueue>
{
public:
  void Init(const TyArrVals &rate_vec, double t0 = 0)
  {
    // each poission_queue_vec[] should have at least one event
    resize(rate_vec.size());
    for (size_t j = 0; j < rate_vec.size(); j++) {
      operator[](j).Init(rate_vec[j], t0);
    }
  }
  TyPoissonTimeVec(const TyArrVals &rate_vec, double t0 = 0)
  {
    Init(rate_vec, t0);
  }
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

  // evolve assume not reset, not external input at all (isolated)
  void NextDtContinuousVec(struct TyNeuronalDymState &nds, double dt)
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

  struct TyNeuronalDymState tmp_neu_state;

  void NextDt(double dt)
  {
    tmp_neu_state = neu_state;
    // evolve as if not reset not resting
    NextDtContinuousVec(tmp_neu_state, dt);
    // keep fired neuron in resting state;
    for (int j = 0; j < pm.n_total(); j++) {
      /// TODO: not exactly correct
      if (tmp_neu_state.time_in_refractory[j]) {
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    TySpikeEventQueue local_spike_events;
    // do reset exceed threshold
    for (int j = 0; j < pm.n_total(); j++) {
      /// TODO: not exactly correct
      if (tmp_neu_state.dym_vals(j, 0) > LIF_param.Vot_Threshold) {
        local_spike_events.emplace(t, j);
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    // do synaptic coupling
    while (local_spike_events.size()) {
      const TySpikeEvent &se = local_spike_events.top();
      if (se.id < pm.n_E) {
        // Excitatory neuron fired
        for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
          if (it.row() < pm.n_E) {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scee;
          } else {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scie;
          }
        }
      } else {
        // Inhibitory neuron fired
        for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
          if (it.row() < pm.n_E) {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scei;
          } else {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scii;
          }
        }
      }
      local_spike_events.pop();
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
      neu_state.dym_vals(se.id, LIF_param.id_gE) += pm.arr_ps[se.id];
//      cout << "Poisson input: to neuron " << se.id
//           << " at " << se.time << endl;
//      cout << "t = " << t << endl;
//      cout << neu_state.dym_vals << endl;
      poisson_events.PopEvent(pm.arr_pr);
    }
    if (t < t_step_end) {
      NextDt(t_step_end - t);
    }
//    cout << "t = " << t << endl;
//    cout << neu_state.dym_vals << endl;

    static int out_i = 0;
    if ((out_i++ & 3)==0) {
      const Eigen::ArrayXXd stat = neu_state.dym_vals.transpose();
      for (int n = 0; n < pm.n_total(); n++) {
        fout << stat.data()+n << '\t';
      }
      fout << '\n';
//      fout << neu_state.dym_vals << endl;
    }
  }
};

class CNeuronSimulatorVer1
{
public:
  struct TyNeuronalParams pm;
  struct TyNeuronalDymState neu_state;
  TyPoissonTimeVec poisson_time_vec;
  double t, dt;

  CNeuronSimulatorVer1(const TyNeuronalParams &_pm, double _dt)
  :pm(_pm), neu_state(_pm)
  {
    t = 0;
    dt = _dt;
    poisson_time_vec.Init(pm.arr_pr, t);
  }

  // Get instantaneous dv/dt for LIF model
  double LIFGetDv(double *dym_val)
  {
    dv_vec = - LIF_param.Con_Leakage * (dym_val(LIF_param.id_V) - LIF_param.Vot_Leakage)
             - dym_val(LIF_param.id_gE) * (dym_val(LIF_param.id_V) - LIF_param.Vot_Excitatory)
             - dym_val(LIF_param.id_gI) * (dym_val(LIF_param.id_V) - LIF_param.Vot_Inhibitory);
  }

  // evolve assume not reset, not external input at all (isolated)
  void NextDtContinuous(double *dym_val, double dt)
  {
    // classical Runge Kutta 4 (for voltage)
    // conductance evolve exactly use exp()

    double tmp_dym_val = dym_val;
    double k1, k2, k3, k4;
    double exp_E = exp(-0.5 * dt / LIF_param.Time_ExCon);
    double exp_I = exp(-0.5 * dt / LIF_param.Time_InCon);

    // k1 = f(t_n, y_n)
    k1 = LIFGetDv(tmp_dym_val);

    // y_n + 0.5*k1*dt
    tmp_dym_val(LIF_param.id_V ) = dym_val(LIF_param.id_V ) + 0.5 * dt * k1;
    tmp_dym_val(LIF_param.id_gE) *= exp_E;
    tmp_dym_val(LIF_param.id_gI) *= exp_I;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k2 = LIFGetDv(tmp_dym_val);

    // y_n + 0.5*k1*dt
    tmp_dym_val(LIF_param.id_V ) = dym_val(LIF_param.id_V ) + 0.5 * dt * k2;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k3 = LIFGetDv(tmp_dym_val);

    // y_n + 0.5*k1*dt
    tmp_dym_val(LIF_param.id_V ) = dym_val(LIF_param.id_V ) + dt * k3;
    tmp_dym_val(LIF_param.id_gE) *= exp_E;
    tmp_dym_val(LIF_param.id_gI) *= exp_I;
    // k2 = f(t+dt/2, y_n + 0.5*k1*dt)
    k4 = LIFGetDv(tmp_dym_val);

    dym_val(LIF_param.id_V ) += dt/6.0 * (k1 + 2*k2 + 2*k3 + k4);
    dym_val(LIF_param.id_gE) = tmp_dym_val(LIF_param.id_gE);
    dym_val(LIF_param.id_gI) = tmp_dym_val(LIF_param.id_gI);
  }

  struct TyNeuronalDymState tmp_neu_state;

  void NextDt(double dt)
  {
    tmp_neu_state = neu_state;
    // evolve as if not reset not resting
    NextDtContinuousVec(tmp_neu_state, dt);
    // keep fired neuron in resting state;
    for (int j = 0; j < pm.n_total(); j++) {
      /// TODO: not exactly correct
      if (tmp_neu_state.time_in_refractory[j]) {
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    TySpikeEventQueue local_spike_events;
    // do reset exceed threshold
    for (int j = 0; j < pm.n_total(); j++) {
      /// TODO: not exactly correct
      if (tmp_neu_state.dym_vals(j, 0) > LIF_param.Vot_Threshold) {
        local_spike_events.emplace(t, j);
        tmp_neu_state.dym_vals(j, 0) = LIF_param.VOT_RESET;
        tmp_neu_state.time_in_refractory[j] += dt;
        if (tmp_neu_state.time_in_refractory[j] > LIF_param.TIME_REFRACTORY) {
          tmp_neu_state.time_in_refractory[j] = 0;
        }
      }
    }
    // do synaptic coupling
    while (local_spike_events.size()) {
      const TySpikeEvent &se = local_spike_events.top();
      if (se.id < pm.n_E) {
        // Excitatory neuron fired
        for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
          if (it.row() < pm.n_E) {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scee;
          } else {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gE) += pm.scie;
          }
        }
      } else {
        // Inhibitory neuron fired
        for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
          if (it.row() < pm.n_E) {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scei;
          } else {
            tmp_neu_state.dym_vals(it.row(), LIF_param.id_gI) += pm.scii;
          }
        }
      }
      local_spike_events.pop();
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
      neu_state.dym_vals(se.id, LIF_param.id_gE) += pm.arr_ps[se.id];
//      cout << "Poisson input: to neuron " << se.id
//           << " at " << se.time << endl;
//      cout << "t = " << t << endl;
//      cout << neu_state.dym_vals << endl;
      poisson_events.PopEvent(pm.arr_pr);
    }
    if (t < t_step_end) {
      NextDt(t_step_end - t);
    }
//    cout << "t = " << t << endl;
//    cout << neu_state.dym_vals << endl;

    static int out_i = 0;
    if ((out_i++ & 3)==0) {
      const Eigen::ArrayXXd stat = neu_state.dym_vals.transpose();
      for (int n = 0; n < pm.n_total(); n++) {
        fout << stat.data()+n << '\t';
      }
      fout << '\n';
//      fout << neu_state.dym_vals << endl;
    }
  }
};

int main()
{
  // fill in parameters
  int n_neu = 100;
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

  for (int i = 0; i < (int)(1e3 / e_dt); i++) {
    neu_simu.NextStep();
  }

  // Octave
  // load a.txt
  // plot([a(1:2:end, 1), a(2:2:end, 1)], '-o')

  return 0;
}
