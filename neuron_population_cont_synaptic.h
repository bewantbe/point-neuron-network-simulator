#include "neuron_population.h"

void hermit(double a, double b, double va, double vb,
            double dva, double dvb, double x, double &fx)
{
  // use v(a), v(b), dv(a)/dt, dv(b)/dt to construct a cubic polynomial
  // basic function f1 satisfies: f1(a)=1, f1(b)=0, f1'(a)=0, f1'(b)=0
  double f1 = va*(2*x+b-3*a)*(x-b)*(x-b)/(b-a)/(b-a)/(b-a);
  // basic function f2 satisfies: f2(a)=0, f2(b)=1, f2'(a)=0, f2'(b)=0
  double f2 = vb*(3*b-2*x-a)*(x-a)*(x-a)/(b-a)/(b-a)/(b-a);
  // basic function f3 satisfies: f3(a)=0, f3(b)=0, f3'(a)=1, f3'(b)=0
  double f3 = dva*(x-a)*(x-b)*(x-b)/(b-a)/(b-a);
  // basic function f4 satisfies: f4(a)=0, f4(b)=0, f4'(a)=0, f4'(b)=1
  double f4 = dvb*(x-a)*(x-a)*(x-b)/(b-a)/(b-a);

  // the polynomial of v(x) - Vot_Threshold
  double V_threshold = 7.5;
  fx = f1 + f2 + f3 + f4 - V_threshold;
}

double root_search(void (*func)(double a, double b, double va, double vb,
                                double dva, double dvb, double x, double &fx),
                   double x1, double x2,
                   double fx1, double fx2,
                   double dfx1, double dfx2, double xacc)
{
//printf("t1, t2 = %.16e, %.16e; v1, v2 = %.16e, %.16e; dv1, dv2 = %.16e, %.16e\n", x1, x2, fx1, fx2, dfx1, dfx2);
  int j;
  double tempx1,tempx2,dx,f,fmid,xmid,root;
  double x1_0 = x1;
  x2 = x2 - x1;
  x1 = 0;

  // for firing time case, fx1<0, fmid>0
  (*func)(x1,x2,fx1,fx2,dfx1,dfx2,x1,f);
  (*func)(x1,x2,fx1,fx2,dfx1,dfx2,x2,fmid);
  /**************************************
    if (fabs(x2-x1)<xacc)
    {
      return x1;
    }
  ***************************************/
  if (f*fmid > 0) {
    fprintf(stderr, "No spike in this interval: v1=%g, v2=%g\n", f, fmid);
    return x1 + x1_0;
  }

  tempx1 = x1;
  tempx2 = x2;
  int Maxnum_search = 50;
  for (j=0; j<Maxnum_search; j++) {
    dx = tempx2 - tempx1;
    xmid = tempx1 + dx/2;
    (*func)(x1,x2,fx1,fx2,dfx1,dfx2,xmid,fmid);
    if (fmid <= 0.0)
      tempx1 = xmid;
    else
      tempx2 = xmid;
    // the interval is small enough or already find the root
    if (fabs(fmid)<xacc) {
      root = xmid;
//  printf("N iter = %d, fmid = %.3e, dx = %.3e, root = %.16e\n", j, fmid, dx, root);
      return root + x1_0;
    }
  }
  return xmid + x1_0;
}


class NeuronPopulationContSyn
  : public NeuronPopulationBase,
    public TyNeuronalParams,
    public TyNeuronalDymState
{
public:
  // Neuron model
  double V_Na = 11.5;
  double V_K  = -1.2;
  double V_L  =  1.06;
  double G_Na = 120;
  double G_K  =  36;
  double G_L  =   0.3;
  double V_gE =  6.5;
  double V_gI = -1.5;
  double tau_gE    = 0.5;
  double tau_gE_s1 = 3.0;
  double tau_gI    = 0.5;
  double tau_gI_s1 = 7.0;
  double V_threshold = 7.5;
  static const int n_var = 8;
  static const int n_var_soma = 4;  // number of variables for non- G part
  static const int id_V     = 0;    // id_V should just before gating variables(see main.cpp)
  static const int id_h     = 6;
  static const int id_m     = 5;
  static const int id_n     = 7;
  static const int id_gE    = 1;
  static const int id_gI    = 3;
  static const int id_gE_s1 = 2;
  static const int id_gI_s1 = 4;
  static const int id_gEInject = id_gE_s1;
  static const int id_gIInject = id_gI_s1;

  void get_dx_net(TyNeuronalDymState::TyDymVals &dx, const TyNeuronalDymState::TyDymVals &xm)
  {
    double synaptic_volt;
    for (int i = 0; i < n_E; i++) {
      if(xm(i, id_V) <= 4) continue;
      synaptic_volt = 1/(1+exp(-5*(xm(i, id_V) - 8.5)));
      for (SparseMat::InnerIterator it(net, i); it; ++it) {
        if (it.row() < n_E) {
          dx(it.row(), id_gEInject) += scee * it.value() * synaptic_volt;
        } else {
          dx(it.row(), id_gEInject) += scie * it.value() * synaptic_volt;
        }
      }
    }
    for (int i = n_E; i < n_neurons(); i++) {
      if(xm(i, id_V) <= 4) continue;
      synaptic_volt = 1/(1+exp(-5*(xm(i, id_V) - 8.5)));
      for (SparseMat::InnerIterator it(net, i); it; ++it) {
        if (it.row() < n_E) {
          dx(it.row(), id_gIInject) += scei * it.value() * synaptic_volt;
        } else {
          dx(it.row(), id_gIInject) += scii * it.value() * synaptic_volt;
        }
      }
    }
  }

  void get_dx(TyNeuronalDymState::TyDymVals &dx,
        const TyNeuronalDymState::TyDymVals &xm, double t)
  {
    const Eigen::ArrayXd &v  = xm.col(id_V);
    const Eigen::ArrayXd &gE = xm.col(id_gE);
    const Eigen::ArrayXd &hE = xm.col(id_gE_s1);
    const Eigen::ArrayXd &gI = xm.col(id_gI);
    const Eigen::ArrayXd &hI = xm.col(id_gI_s1);
    const Eigen::ArrayXd &m  = xm.col(id_m);
    const Eigen::ArrayXd &h  = xm.col(id_h);
    const Eigen::ArrayXd &n  = xm.col(id_n);

    // cost 13.3% (1000neuron), 16.4% (100neuron)
    dx.col(id_V) = -G_L * (v - V_L)
                -G_Na * m*m*m*h * (v-V_Na)
                -G_K * n*n*n*n * (v-V_K)
                -gE * (v-V_gE) - gI * (v-V_gI);

    // cost 55.6% (1000neuron), 59.8% (100neuron) - 52.0% (1000neuron), % (100neuron)
    // mhn
    static Eigen::ArrayXd e1, e2;
    if (e1.size() != n_neurons()) {
      e1.resize(n_neurons());
      e2.resize(n_neurons());
    }
    e1 = 2.5-v;
    e1 = exp(e1);     // only write in this way, Eigen will use MKL
    e2 = v*(-1/1.8);
    e2 = exp(e2);
    dx.col(5) = (v-2.5)/(1-e1)*(1-m) - 4*e2*m;
    e2 = e1*exp(-2.5);
    e2 = sqrt(e2);
    dx.col(6) = 0.07*e2*(1-h) - h/(1+e1*exp(0.5));
    e2 = sqrt(e2);
    e2 = sqrt(e2);
    dx.col(7) = 0.1*(v-1.0)/(1-e1*exp(-1.5))*(1-n) - 0.125*e2*n;

    //conductance_dt
    dx.col(1) = -gE * (1.0/tau_gE) + hE;
    dx.col(2) = -hE * (1.0/tau_gE_s1);
    dx.col(3) = -gI * (1.0/tau_gI) + hI;
    dx.col(4) = -hI * (1.0/tau_gI_s1);

    // cost 12% (1000neuron), 2.4% (100neuron)
    // Compute influence from the network (H^E and H^I)
    get_dx_net(dx, xm);
  }

  void voltage_dt(int index_neuron,double t,double m,double h,double n,double gE,double gI,double v,double &dv)
  {
    dv = - G_L*(v - V_L) - G_Na*m*m*m*h*(v-V_Na) - G_K*n*n*n*n*(v - V_K)
            - gE*(v - V_gE) - gI*(v - V_gI);
  }

  void spiking_time(double *bef_neu, double *cur_neu, double t_evolution, double Ta,
                    double Tb, double &firing_time)
  {
    double dva = 0; // initialization
    double va  = bef_neu[id_V];
    double gEa = bef_neu[id_gE];
    double gIa = bef_neu[id_gI];
    double ma  = bef_neu[id_m];
    double ha  = bef_neu[id_h];
    double na  = bef_neu[id_n];

    voltage_dt(0,t_evolution,ma,ha,na,gEa,gIa,va,dva);

    double dvb = 0; // initialization
    double vb  = cur_neu[0];
    double gEb = cur_neu[1];
    double gIb = cur_neu[3];
    double mb  = cur_neu[5];
    double hb  = cur_neu[6];
    double nb  = cur_neu[7];

    voltage_dt(0,(t_evolution+Tb-Ta),mb,hb,nb,gEb,gIb,vb,dvb);

// the tolerance error for root searching
#define root_acc  (1.0e-12)

    firing_time = root_search(hermit, Ta, Tb, va, vb, dva, dvb, root_acc);
  }

  TyNeuronalDymState::TyDymVals xt, k1, k2, k3, k4;  // all variables of all neurons each
  void runge_kutta4_vec(TyNeuronalDymState::TyDymVals &neu_dym_vals,
      double dt, double t_local, TySpikeEventVec &spike_events)
  {
    const int n = neu_dym_vals.rows();
    const int m = neu_dym_vals.cols();
    if (xt.rows() != n || xt.cols() != m) {
      xt.resize(n, m);
      k1.resize(n, m);
      k2.resize(n, m);
      k3.resize(n, m);
      k4.resize(n, m);
    }

    xt = neu_dym_vals;

    get_dx(k1, xt, t_local);

    get_dx(k2, xt + k1 * (dt/2), t_local + dt/2);

    get_dx(k3, xt + k2 * (dt/2), t_local + dt/2);

    get_dx(k4, xt + k3 * dt, t_local + dt);

    neu_dym_vals = xt + (k1 + 2*k2 + 2*k3 + k4) * (dt/6);
    
    // Get spike timings
    size_t id_ses_hd = spike_events.size();
    for (int j=0; j<n; j++) {
      if (   neu_dym_vals(j, id_V) >= V_threshold
                    && xt(j, id_V) <  V_threshold) {
        double t_sp = NAN;
        spiking_time(xt.data()+j*m, neu_dym_vals.data()+j*m,
                     t_local, t_local, t_local+dt, t_sp);
        spike_events.emplace_back(t_sp, j);
      }
    }
    std::sort(spike_events.begin() + id_ses_hd, spike_events.end());
  }

  // ---------------------------------------------

  NeuronPopulationContSyn(const TyNeuronalParams &_pm)
    :TyNeuronalParams(_pm), TyNeuronalDymState(_pm, n_var)
  {
  }

  int n_neurons() const
  {
    return Get_n_neurons();
  }

  void NoInteractDt(int neuron_id, double dt, double t_local,
                    TySpikeEventVec &spike_events)
  {
    fprintf(stderr, "Evoluate each neuron separately is not supported.\n");
    throw "Evoluate each neuron separately is not supported.";
  }

  void NoInteractDt(double dt, double t_local,
                    TySpikeEventVec &spike_events)
  {
    runge_kutta4_vec(dym_vals, dt, t_local, spike_events);
  }

  void ForceReset(int neuron_id)
  {
  }

  void InjectPoissonE(int neuron_id)
  {
    // Add Poisson input
    dym_vals(neuron_id, id_gEInject) += arr_ps[neuron_id];
  }

  void operator=(const TyNeuronalDymState &neu_dym)
  {
    GetDymState() = neu_dym;
  }

  void ScatterCopy(const struct TyNeuronalDymState &nd,
                   const std::vector<int> &ids)
  {
    GetDymState().ScatterCopy(nd, ids);
  }

  TyNeuronalDymState & GetDymState()
  {
    return *(static_cast<TyNeuronalDymState*>(this));
  }
  const TyNeuronalDymState & GetDymState() const
  {
    return *(static_cast<const TyNeuronalDymState*>(this));
  }
  TyNeuronalDymState * GetDymStatePtr()
  {
    return static_cast<TyNeuronalDymState*>(this);
  }
  const TyNeuronalDymState * GetDymStatePtr() const
  {
    return static_cast<const TyNeuronalDymState*>(this);
  }
  double GetDymState(int neuron_id, int id_dym) const
  {
    return dym_vals(neuron_id, id_dym);
  }
  const TyNeuronalParams * GetNeuronalParamsPtr() const
  {
    return static_cast<const TyNeuronalParams*>(this);
  }

  void SynapticInteraction(const TySpikeEvent &se) override
  {}
  void SynapticInteraction(int neuron_id, const TySpikeEvent &se) override
  {}
};
