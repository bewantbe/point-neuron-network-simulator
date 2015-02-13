#include "common_header.h"
#include "math_helper.h"

// #define NDEBUG
#include <cassert>

const bool g_b_debug = false;
//#define dbg_printf(...) ((void)0);
#define dbg_printf printf

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

  TyPoissonTimeVec() = default;

  void RestoreIdx()
  {
    for (size_t j = 0; j < size(); j++) {
      operator[](j).id_seq = id_seq_vec[j];
    }
  }

  void SaveIdxAndClean()
  {
    int j = 0;
    for (iterator it = begin(); it != end(); it++, j++) {
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

    if (v_n <= LIF_param.Vot_Threshold
        && dym_val[LIF_param.id_V ] > LIF_param.Vot_Threshold) {
      spike_time_local = root_search(0, dt,
        v_n, dym_val[LIF_param.id_V ],
        k1, LIFGetDv(dym_val), LIF_param.Vot_Threshold, dt, 1e-12);
    } else {
      if (v_n > 0.996 && k1>0) { // the v_n > 0.996 is for dt=0.5 ms
        // try extract some missing spikes
        double c = v_n;
        double b = k1;
        double a = (dym_val[LIF_param.id_V ] - c - b*dt)/(dt*dt);
        double t_max_guess = -b/(2*a);
        // in LIF, it can guarantee that a<0 (concave), hence t_max_guess > 0
        if (t_max_guess<dt) {
          dbg_printf("rare case captured, guess time: %f\n", t_max_guess);
          spike_time_local = root_search(0, dt,
            v_n, dym_val[LIF_param.id_V ],
            k1, LIFGetDv(dym_val), LIF_param.Vot_Threshold, t_max_guess, 1e-12);
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
    dym_val[LIF_param.id_gE] *= exp(-dt / LIF_param.Time_ExCon);
    dym_val[LIF_param.id_gI] *= exp(-dt / LIF_param.Time_InCon);
  }

  void SynapticInteraction(struct TyNeuronalDymState &e_neu_state, const TySpikeEvent &se) const
  {
    if (se.id < pm.n_E) {
      // Excitatory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gE) += it.value() * pm.scee;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gE) += it.value() * pm.scie;
        }
      }
    } else {
      // Inhibitory neuron fired
      for (SparseMat::InnerIterator it(pm.net, se.id); it; ++it) {
        if (it.row() < pm.n_E) {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gI) += it.value() * pm.scei;
        } else {
          e_neu_state.dym_vals(it.row(), LIF_param.id_gI) += it.value() * pm.scii;
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
        // this std::max is used to avoid very rare case
        t_in_refractory = std::max(dt_local - spike_time_local,
                                   std::numeric_limits<double>::min());
        if (t_in_refractory < LIF_param.TIME_REFRACTORY) {
          dym_val[LIF_param.id_V] = LIF_param.VOT_RESET;
        } else {
          // short refractory period (< dt_local), active again
          dt_local = t_in_refractory - LIF_param.TIME_REFRACTORY;
          t_in_refractory = 0;
          // set conductance state at firing
          NextDtConductance(dym_val, -dt_local);
          dym_val[LIF_param.id_V] = LIF_param.VOT_RESET;
          double spike_time_local_tmp;
          NextDtContinuous(dym_val, spike_time_local_tmp, dt_local);
          if (!std::isnan(spike_time_local_tmp)) {
            cerr << "NextQuietDtNeuState(): spike too fast" << endl;
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
          cerr << "NextQuietDtNeuState(): spike too fast" << endl;
        }
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
        NextDtConductance(dym_val, dt_local);
        t_in_refractory += dt_local;
      }
    }
  }

  // evolve as if no synaptic interaction
  void NextDt(struct TyNeuronalDymState &tmp_neu_state,
              TySpikeEventQueue &spike_events, double dt)
  {
    double t_step_end = t + dt;
    dbg_printf("----- Trying t = %f .. %f -------\n", t, t_step_end);
    // evolve each neuron as if no reset
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
    double t_end = t + dt;
    TySpikeEvent latest_spike_event =
      {std::numeric_limits<double>::quiet_NaN(), -1};
    struct TyNeuronalDymState bk_neu_state;

    dbg_printf("===== NextStep(): t = %f .. %f\n", t, t_end);
    poisson_time_vec.SaveIdxAndClean();
    do {
      TySpikeEventQueue spike_events;
      bk_neu_state = neu_state;         // make a backup
      //poisson_time_vec.RestoreIdx();
      NextDt(neu_state, spike_events, t_end - t);
      if (spike_events.size() != 0) {
        neu_state = bk_neu_state;
        dbg_printf("Neuron [%d] spike at t = %f\n",
               spike_events.top().id, spike_events.top().time);
        latest_spike_event = spike_events.top();
        // Really evolve the whole system. spike_events will be discard.
        poisson_time_vec.RestoreIdx();
        NextDt(neu_state, spike_events, latest_spike_event.time - t);
        if (spike_events.top().id != latest_spike_event.id) {
          dbg_printf("error: NextStep(): spike neuron is different from what it should be!\n  Original: #%d at %f, now #%d at %f\n",
              latest_spike_event.id, latest_spike_event.time,
              spike_events.top().id, spike_events.top().time) ;
          dbg_printf("spike events in this interval:\n");
          while (spike_events.size() > 0) {
            TySpikeEvent se = spike_events.top();
            dbg_printf("  #%2d at %f, v0=%f g0=%f, v1=%f g1=%f\n",
            se.id, se.time,
            bk_neu_state.dym_vals(se.id, 0), bk_neu_state.dym_vals(se.id, 1),
            neu_state.dym_vals(se.id, 0), neu_state.dym_vals(se.id, 1));
            spike_events.pop();
          }
        }
        t = latest_spike_event.time;
        // force the neuron to spike, if not already
        if (neu_state.time_in_refractory[latest_spike_event.id] == 0) {
          neu_state.dym_vals(latest_spike_event.id, LIF_param.id_V) = LIF_param.VOT_RESET;
          neu_state.time_in_refractory[latest_spike_event.id] =
            std::numeric_limits<double>::min();
        }
        SynapticInteraction(neu_state, latest_spike_event);
        poisson_time_vec.SaveIdxAndClean();
      } else {
        break;
      }
    } while (true);
    t = t_end;
  }
};

typedef CNeuronSimulatorEvolveEach CNeuronSimulator;

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

int main(int argc, char *argv[])
{
  // fill in parameters
  int n_neu = 2;
  TyNeuronalParams pm(LIF_G, n_neu, 0);

  cout << "network:\n";
  cout << pm.net << "\n";
  cout << "end network\n";

  for (int i = 0; i < n_neu; i++) {
    pm.arr_pr[i] = 1.0;
    pm.arr_ps[i] = 0.012;
  }
  typedef Eigen::Triplet<double> TyEdgeTriplet;
  std::vector< TyEdgeTriplet > net_coef;
  //net_coef.push_back({1, 0, 1.0});
  for (int i = 0; i < n_neu; i++) {
    for (int j = 0; j < n_neu; j++) {
      if (i==j) {
        continue;
      }
      net_coef.push_back(TyEdgeTriplet(i, j, 1.0));
    }
  }
  cout << "size net_coef = " << net_coef.size() << "\n";
  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());
  for (int j = 0; j < pm.n_total(); j++) {
    if (pm.net.coeffRef(j,j)) {
      pm.net.coeffRef(j,j) = 0;
    }
  }
  pm.net.prune(std::numeric_limits<double>::min(), 1);  // remove zeros;
  pm.net.makeCompressed();

  for (int k=0; k<2; k++) {
    for (SparseMat::InnerIterator it(pm.net, k); it; ++it) {
      cout << "edge (" << it.row() << " " << k << ")" << endl;
    }
  }

  cout << "network:\n";
  cout << pm.net << "\n";
  cout << "end network\n";

  //return 1;

  pm.scee = 0.01;
  pm.scie = 0;
  pm.scei = 0;
  pm.scii = 0;

  // Create simulator
  double e_dt = 1.0 / 2;  // ms
  CNeuronSimulator neu_simu(pm, e_dt);

  FillNeuStateFromFile(neu_simu.neu_state, "neu_state_init.txt");
  FillPoissonEventsFromFile(neu_simu.poisson_time_vec, "poisson_events.txt");
  for (int i = 0; i < pm.n_total(); i++) {
    cout << "neuron [" << i << "] poisson vec " << neu_simu.poisson_time_vec[i].size() << endl;
  }

  cout << "t = " << neu_simu.t << endl;
  cout << neu_simu.neu_state.dym_vals << endl;

  FILE *volt_fout = fopen("volt.txt", "w");

  for (int j = 0; j < pm.n_total(); j++) {
    fprintf(volt_fout, "%.16e %.16e %.16e  ",
      neu_simu.neu_state.dym_vals(j, 0),
      neu_simu.neu_state.dym_vals(j, 1),
      neu_simu.neu_state.dym_vals(j, 2));
  }
  fprintf(volt_fout, "\n");

  //std::ofstream volt_fout("volt.txt");
  for (int i = 0; i < (int)(1e3 / e_dt); i++) {
    neu_simu.NextStep();
    for (int j = 0; j < pm.n_total(); j++) {
      fprintf(volt_fout, "%.16e %.16e %.16e  ",
        neu_simu.neu_state.dym_vals(j, 0),
        neu_simu.neu_state.dym_vals(j, 1),
        neu_simu.neu_state.dym_vals(j, 2));
    }
    fprintf(volt_fout, "\n");
    if (g_b_debug) {
      cout << "t = " << neu_simu.t << endl;
      cout << neu_simu.neu_state.dym_vals << endl;
    }
  }
  fclose(volt_fout);

  cout << "t = " << neu_simu.t << endl;
  cout << neu_simu.neu_state.dym_vals << endl;

  // Octave
  // load a.txt
  // plot([a(1:2:end, 1), a(2:2:end, 1)], '-o')

  return 0;
}
