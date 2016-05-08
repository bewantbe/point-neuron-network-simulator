#ifndef HEADER_SINGLE_NEURON_DYNAMICS
#define HEADER_SINGLE_NEURON_DYNAMICS

#include "common_header.h"
#include "math_helper.h"

/*
// Fast code for exp(x)
#include "fmath.hpp"
//#define exp(x) fmath::expd(x)

template<typename Ty>
inline Ty my_expd(const Ty &x)
{ return exp(x); }

template<>
inline double my_expd(const double &x)
{ return fmath::expd(x); }

#define exp(x) my_expd(x)
*/

struct Ty_Neuron_Dym_Base
{
  virtual double Get_V_threshold() const = 0;
  virtual int Get_id_V() const = 0;
  virtual int Get_id_gE() const = 0;
  virtual int Get_id_gI() const = 0;
  virtual int Get_id_gEInject() const = 0;  // index of gE injection variable
  virtual int Get_id_gIInject() const = 0;  // index of gI injection variable
  virtual int Get_n_dym_vars() const = 0;   // number of state variables
  virtual void NextStepSingleNeuronQuiet(
    double *dym_val,
    double &t_in_refractory,
    double &spike_time_local,
    double dt_local) const = 0;
  virtual void VoltHandReset(double *dym_val) const = 0;
};

struct Ty_LIF_G_core
{
  // The neuron model named LIF-G in this code.
  // This is the Leaky Integrate-and-Fire model with jump conductance.
  double V_threshold   = 1.0;  // voltages are in dimensionless unit
  double V_reset       = 0.0;
  double V_leakage     = 0.0;
  double V_excitatory  = 14.0/3.0;
  double V_inhibitory  = -2.0/3.0;
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
    return - Con_Leakage * (dym_val[id_V] - V_leakage)
           - dym_val[id_gE] * (dym_val[id_V] - V_excitatory)
           - dym_val[id_gI] * (dym_val[id_V] - V_inhibitory);
  }

  // Evolve the state `dym_val' a `dt' forward,
  // using classical Runge–Kutta 4-th order scheme for voltage.
  // Conductance will evolve using the exact formula.
  // Return derivative k1 at t_n, for later interpolation.
  MACRO_NO_INLINE double DymInplaceRK4(double *dym_val, double dt) const
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

struct Ty_LIF_GH_core
{
  // The neuron model named LIF-GH in this code.
  // This is the Leaky Integrate-and-Fire model with order 1 smoothed conductance.
  double V_threshold   = 1.0;  // voltages are in dimensionless unit
  double V_reset       = 0.0;
  double V_leakage     = 0.0;
  double V_excitatory  = 14.0/3.0;
  double V_inhibitory  = -2.0/3.0;
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
  inline double GetDv(const double *dym_val) const
  {
    return - Con_Leakage * (dym_val[id_V] - V_leakage)
           - dym_val[id_gE] * (dym_val[id_V] - V_excitatory)
           - dym_val[id_gI] * (dym_val[id_V] - V_inhibitory);
  }

  // Evolve the state `dym_val' a `dt' forward,
  // using classical Runge–Kutta 4-th order scheme for voltage.
  // Conductance will evolve using the exact formula.
  // Return derivative k1 at t_n, for later interpolation.
  MACRO_NO_INLINE double DymInplaceRK4(double *dym_val, double dt) const
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

// Adapter for IF type models (only sub-threshold dynamics and need hand reset)
template<typename TyNeuronModel>
struct Ty_LIF_stepper: public TyNeuronModel, public Ty_Neuron_Dym_Base
{
  // for template class we need these "using"s. It's a requirement for TyNeuronModel.
  using TyNeuronModel::id_V;
  using TyNeuronModel::id_gE;
  using TyNeuronModel::V_threshold;
  using TyNeuronModel::V_reset;
  using TyNeuronModel::Time_Refractory;
  using TyNeuronModel::GetDv;
  using TyNeuronModel::NextDtConductance;
  using TyNeuronModel::DymInplaceRK4;

  using TyNeuronModel::id_gEInject;
  using TyNeuronModel::id_gIInject;
  using TyNeuronModel::n_var;
  using TyNeuronModel::id_gI;
  double Get_V_threshold() const override {return V_threshold;};
  int Get_id_gEInject() const override {return id_gEInject;}
  int Get_id_gIInject() const override {return id_gIInject;}
  int Get_n_dym_vars() const override {return n_var;}
  int Get_id_V() const override {return id_V;}
  int Get_id_gE() const override {return id_gE;}
  int Get_id_gI() const override {return id_gI;}

  // Used when reset the voltage by hand. (e.g. outside this class)
  inline void VoltHandReset(double *dym_val) const override
  {
    dym_val[id_V] = V_reset;
  }

  // Evolve the ODE and note down the spike time, assuming no reset and no external input.
  // `spike_time_local' should be guaranteed to be within [0, dt] or NAN.
  MACRO_NO_INLINE void NextStepSingleNeuronContinuous(double *dym_val, double &spike_time_local, double dt) const
  {
    double v_n = dym_val[id_V];
    double k1  = DymInplaceRK4(dym_val, dt);

    if (v_n <= V_threshold
        && dym_val[id_V ] > V_threshold) {
      spike_time_local = cubic_hermit_real_root(dt,
        v_n, dym_val[id_V ],
        k1, GetDv(dym_val), V_threshold);
    } else {
      if (v_n > 0.996 && k1>0) { // the v_n > 0.996 is for dt=0.5 ms, LIF,G model
        // Try capture some missing spikes that the intermediate value passes
        // threshold, but both ends are lower than threshold.
        // Get quadratic curve from value of both ends and derivative from left end
        // Return the maximum point as `t_max_guess'
        double c = v_n;
        double b = k1;
        double a = (dym_val[id_V ] - c - b*dt)/(dt*dt);
        double t_max_guess = -b/(2*a);
        // In LIF-G, it can guarantee that a<0 (concave),
        // hence t_max_guess > 0. But in LIF-G model, we still need to
        // check 0 < t_max_guess
        if (0 < t_max_guess && t_max_guess < dt
            && (b*b)/(-4*a)+c >= V_threshold) {
          //dbg_printf("Rare event: mid-dt spike captured, guess time: %f\n", t_max_guess);
          dbg_printf("NextStepSingleNeuronContinuous(): possible mid-dt spike detected:\n");
          dbg_printf("  Guessed max time: %f, dt = %f\n", t_max_guess, dt);
          // root should in [0, t_max_guess]
          spike_time_local = cubic_hermit_real_root(dt,
            v_n, dym_val[id_V ],
            k1, GetDv(dym_val), V_threshold);
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
  void NextStepSingleNeuronQuiet(double *dym_val, double &t_in_refractory,
                           double &spike_time_local, double dt_local) const override
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
        dym_val[id_V] = V_reset;
        dbg_printf("Reach threshold detected. Neuron reseted.\n");
        if (t_in_refractory >= Time_Refractory) {
          // Short refractory period (< dt_local), neuron will active again.
          dt_local = t_in_refractory - Time_Refractory;
          t_in_refractory = 0;
          // Back to the activation time.
          NextDtConductance(dym_val, -dt_local);
          double spike_time_local_tmp;
          NextStepSingleNeuronContinuous(dym_val, spike_time_local_tmp, dt_local);
          if (!std::isnan(spike_time_local_tmp)) {
            cerr << "NextStepSingleNeuronQuiet(): Multiple spikes in one dt. Interaction dropped." << endl;
            cerr << "  dt_local = " << dt_local << '\n';
            cerr << "  spike_time_local = " << spike_time_local << '\n';
            cerr << "  t_in_refractory = " << t_in_refractory << '\n';
            cerr << "  dym_val[id_V] = " << dym_val[id_V] << '\n';
            cerr << "  dym_val[id_gE] = " << dym_val[id_gE] << '\n';
            throw "Multiple spikes in one dt.";
          }
        }
      }
    } else {
      // Neuron in refractory period.
      double dt_refractory_remain = Time_Refractory
                                 - t_in_refractory;
      if (dt_refractory_remain < dt_local) {
        // neuron will awake after dt_refractory_remain which is in this dt_local
        NextDtConductance(dym_val, dt_refractory_remain);
        assert( dym_val[id_V] == V_reset );
        t_in_refractory = 0;
        NextStepSingleNeuronQuiet(dym_val, t_in_refractory,
                            spike_time_local, dt_local - dt_refractory_remain);
      } else {
        spike_time_local = std::numeric_limits<double>::quiet_NaN();
        NextDtConductance(dym_val, dt_local);
        t_in_refractory += dt_local;
      }
    }
  }

};

typedef Ty_LIF_stepper<Ty_LIF_G_core>  Ty_LIF_G;
typedef Ty_LIF_stepper<Ty_LIF_GH_core> Ty_LIF_GH;

// Model of classical Hodgkin-Huxley (HH) neuron,
// with two ODE for G (See Ty_LIF_GH_core::NextDtConductance() for details).
// The parameters for G comes from unknown source.
template<typename ExtraCurrent>
struct Ty_HH_GH_CUR
  :public Ty_Neuron_Dym_Base
{
  typedef typename ExtraCurrent::TyData TyCurrentData;

  double V_Na = 115;             // mV
  double V_K  = -12;
  double V_L  =  10.6;
  double G_Na = 120;             // mS cm^-2
  double G_K  =  36;
  double G_L  =   0.3;
  double V_gE =  65;             // mV, for synaptic
  double V_gI = -15;             // mV
  double tau_gE    = 0.5;        // ms, decay time for gE equation
  double tau_gE_s1 = 3.0;        // ms, decay time for gE_H equation.
  double tau_gI    = 0.5;        // ms
  double tau_gI_s1 = 7.0;        // ms
  double V_threshold = 50;       // mV, determine the spike time for synaptic interaction.
  //double Time_Refractory = 1.0;  // ms, hard refractory period. Used for correctly locate the firing when change time step
  static const int n_var = 8;
  static const int n_var_soma = 4;  // number of variables for non- G part
  static const int id_V     = 0;    // id_V should just before gating variables(see main.cpp)
  static const int id_h     = 1;
  static const int id_m     = 2;
  static const int id_n     = 3;
  static const int id_gE    = 4;
  static const int id_gI    = 5;
  static const int id_gE_s1 = 6;
  static const int id_gI_s1 = 7;
  static const int id_gEInject = id_gE_s1;
  static const int id_gIInject = id_gI_s1;

  double Get_V_threshold() const override {return V_threshold;};
  int Get_id_gEInject() const override {return id_gEInject;}
  int Get_id_gIInject() const override {return id_gIInject;}
  int Get_n_dym_vars() const override {return n_var;}
  int Get_id_V() const override {return id_V;}
  int Get_id_gE() const override {return id_gE;}
  int Get_id_gI() const override {return id_gI;}

  // dV/dt term
  inline double GetDv(const double *dym_val, double t,
               const TyCurrentData &extra_data) const
  {
    const double &V = dym_val[id_V];
    const double &h = dym_val[id_h];
    const double &m = dym_val[id_m];
    const double &n = dym_val[id_n];
    /*printf("new. extra_data: %e\n", ExtraCurrent()(t, extra_data));*/
    return
      -(V-V_Na) * G_Na * h * m * m * m
      -(V-V_K ) * G_K * n * n * n * n
      -(V-V_L ) * G_L
      -(V-V_gE) * dym_val[id_gE]
      -(V-V_gI) * dym_val[id_gI]
      + ExtraCurrent()(t, extra_data);
  }

//  // Classical HH neuron equations go here
//  void ODE_RHS(const double *dym_val, double *dym_d_val) const
//  {
//    const double &V = dym_val[id_V];
//    const double &h = dym_val[id_h];
//    const double &m = dym_val[id_m];
//    const double &n = dym_val[id_n];
//    dym_d_val[id_V] = GetDv(dym_val);
//    dym_d_val[id_h] = (1-h) * (0.07*exp(-V/20))
//                        - h / (exp(3-0.1*V)+1);
//    dym_d_val[id_m] = (1-m) * ((2.5-0.1*V)/(exp(2.5-0.1*V)-1))
//                        - m * (4*exp(-V/18));
//    dym_d_val[id_n] = (1-n) * ((0.1-0.01*V)/(exp(1-0.1*V)-1))
//                        - n * (0.125*exp(-V/80));
//    // No RHS for G here, that is solved explicitly.
//  }

  // Classical HH neuron equations go here. Faster version.
  void ODE_RHS(const double *dym_val, double *dym_d_val, double t,
               const TyCurrentData &extra_data) const
  {
    const double &V = dym_val[id_V];
    const double &h = dym_val[id_h];
    const double &m = dym_val[id_m];
    const double &n = dym_val[id_n];
    dym_d_val[id_V] = GetDv(dym_val, t, extra_data);
    double e1 = exp(-0.1*dym_val[id_V]);
    double e2 = sqrt(e1);
    dym_d_val[id_h] = (1-h) * (0.07*e2)
                        - h / (e1*exp(3.0)+1);
    dym_d_val[id_m] = (1-m) * ((2.5-0.1*V)/(e1*exp(2.5)-1))
                        - m * (4*exp(-V/18));
    dym_d_val[id_n] = (1-n) * ((0.1-0.01*V)/(e1*exp(1.0)-1))
                        - n * (0.125*sqrt(sqrt(e2)));
    // No RHS for G here, that is solved explicitly.
  }

  double DymInplaceRK4(double *dym_val, double dt, double t,
                       const TyCurrentData &extra_data) const
  {
    // for G. See Ty_LIF_GH_core::NextDtConductance().
    double expEC  = exp(-0.5 * dt / tau_gE);
    double expECR = exp(-0.5 * dt / tau_gE_s1);
    double expIC  = exp(-0.5 * dt / tau_gI);
    double expICR = exp(-0.5 * dt / tau_gI_s1);
    double gE_s_coef = (expEC - expECR) * tau_gE * tau_gE_s1 / (tau_gE - tau_gE_s1);
    double gI_s_coef = (expIC - expICR) * tau_gI * tau_gI_s1 / (tau_gI - tau_gE_s1);

    double k1[n_var_soma], k2[n_var_soma], k3[n_var_soma], k4[n_var_soma];
    double dym_val0[n_var_soma];
    // Use template BLAS lib, so some expressions can be written in vector form.
    typedef Eigen::Map< Eigen::RowVectorXd > TyMapVec;
    TyMapVec dym_val0_v(dym_val0, n_var_soma);
    TyMapVec dym_val_v (dym_val , n_var_soma);
    TyMapVec k1_v(k1, n_var_soma);
    TyMapVec k2_v(k2, n_var_soma);
    TyMapVec k3_v(k3, n_var_soma);
    TyMapVec k4_v(k4, n_var_soma);

    dym_val0_v = dym_val_v;
    ODE_RHS(dym_val, k1, t, extra_data);

    dym_val_v = dym_val0_v + 0.5*dt*k1_v;
    dym_val[id_gE] = expEC*dym_val[id_gE] + gE_s_coef*dym_val[id_gE_s1];
    dym_val[id_gI] = expIC*dym_val[id_gI] + gI_s_coef*dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    ODE_RHS(dym_val, k2, t + 0.5*dt, extra_data);

    dym_val_v = dym_val0_v + 0.5*dt*k2_v;
    ODE_RHS(dym_val, k3, t + 0.5*dt, extra_data);

    dym_val_v = dym_val0_v + dt*k3_v;
    dym_val[id_gE] = expEC*dym_val[id_gE] + gE_s_coef*dym_val[id_gE_s1];
    dym_val[id_gI] = expIC*dym_val[id_gI] + gI_s_coef*dym_val[id_gI_s1];
    dym_val[id_gE_s1] *= expECR;
    dym_val[id_gI_s1] *= expICR;
    ODE_RHS(dym_val, k4, t + dt, extra_data);

    dym_val_v = dym_val0_v + dt/6.0 * (k1_v + 2*k2_v + 2*k3_v + k4_v);

    return k1[id_V];
  };

  void VoltHandReset(double *dym_val) const override
  {
    // no force reset required for HH
  }

  void NextStepSingleNeuronQuiet(
    double *dym_val,
    double &t_in_refractory,
    double &spike_time_local,
    double dt_local,
    double t,
    const TyCurrentData &extra_data) const
  {
    /*printf("new. extra_data: %e\n", ExtraCurrent()(t, extra_data));*/
    spike_time_local = std::numeric_limits<double>::quiet_NaN();
    double v0 = dym_val[id_V];
    double k1 = DymInplaceRK4(dym_val, dt_local, t, extra_data);
    double &v1 = dym_val[id_V];
    // See if neuron is firing. t_in_refractory == 0 means the neuron
    // is not in hand set refractory period, avoids kind of infinite loop.
    if (v0 <= V_threshold && v1 >= V_threshold && t_in_refractory == 0) {
      spike_time_local = cubic_hermit_real_root(dt_local,
        v0, v1, k1, GetDv(dym_val, t, extra_data), V_threshold);
      // Some redundant calculation of GetDv(). Let's leave it.
    }
    if (dt_local>0) {
      // not 100% mathematically safe. Use the hard refractory for that.
      t_in_refractory = 0;
    }
    //t_in_refractory += dt_local;
    //if (t_in_refractory > Time_Refractory) {
    //t_in_refractory = 0;
    //}
  }
};

// Model that note down spike event when the spike is falling.
template<typename ExtraCurrent>
struct Ty_HH_FT_GH_CUR  // Falling threshold
  :public Ty_HH_GH_CUR<ExtraCurrent>
{
  typedef typename Ty_HH_GH_CUR<ExtraCurrent>::TyCurrentData TyCurrentData;
  using Ty_HH_GH_CUR<ExtraCurrent>::id_V;
  using Ty_HH_GH_CUR<ExtraCurrent>::V_threshold;
  using Ty_HH_GH_CUR<ExtraCurrent>::DymInplaceRK4;
  using Ty_HH_GH_CUR<ExtraCurrent>::GetDv;

  void NextStepSingleNeuronQuiet(
    double *dym_val,
    double &t_in_refractory,
    double &spike_time_local,
    double dt_local,
    double t,
    const TyCurrentData &extra_data) const
  {
    spike_time_local = std::numeric_limits<double>::quiet_NaN();
    double v0 = dym_val[id_V];
    double k1 = DymInplaceRK4(dym_val, dt_local, t, extra_data);
    double &v1 = dym_val[id_V];
    if (v0 >= V_threshold && v1 <= V_threshold && t_in_refractory == 0) {
      spike_time_local = cubic_hermit_real_root(dt_local,
        v0, v1, k1, GetDv(dym_val, t, extra_data), V_threshold);
    }
    if (dt_local>0) {
      t_in_refractory = 0;
    }
  }
};

// Template for zero current input
struct TyZeroCurrent
{
  typedef int TyData;
  double operator()(double t, const TyData &a) const
  { return 0.0; }
};

template<class Ty_Neuron_With_Current>
struct Neuron_Zero_Current_Adaper :public Ty_Neuron_With_Current
{
  /*using Ty_HH_GH_CUR<TyZeroCurrent>::NextStepSingleNeuronQuiet;*/
  void NextStepSingleNeuronQuiet(double *dym_val, double &t_in_refractory,
    double &spike_time_local, double dt_local) const override
  {
    Ty_Neuron_With_Current::NextStepSingleNeuronQuiet(dym_val, t_in_refractory, spike_time_local, dt_local, 0, 0);
  }
};

// Declare the model with zero current input
typedef Neuron_Zero_Current_Adaper< Ty_HH_GH_CUR<TyZeroCurrent> >  Ty_HH_GH;
typedef Neuron_Zero_Current_Adaper< Ty_HH_FT_GH_CUR<TyZeroCurrent> >  Ty_HH_FT_GH;

// Template for sine current input
struct TySineCurrent
{
  typedef double * TyData;  // Amplitude, angular frequency, phase
  double operator()(double t, const TyData &a) const
  { //printf("Current: %e, %e, %e\n", a[0], a[1], a[2]);  fflush(stdout);
    return a[0]*sin(a[1]*t + a[2]); }
};

template<class Ty_Neuron_With_Current>
struct Neuron_Sine_Current_Adaper :public Ty_Neuron_With_Current
{
  using Ty_Neuron_With_Current::NextStepSingleNeuronQuiet;

  void NextStepSingleNeuronQuiet(double *dym_val, double &t_in_refractory,
    double &spike_time_local, double dt_local) const override
  {
    printf("Old.\n");
    double d[3] = {0, 0, 0};
    Ty_Neuron_With_Current::NextStepSingleNeuronQuiet(dym_val, t_in_refractory, spike_time_local, dt_local, 0, d);
  }
};

// Declare the model with sine current input
typedef Neuron_Sine_Current_Adaper< Ty_HH_GH_CUR<TySineCurrent> > Ty_HH_GH_sine;
typedef Neuron_Sine_Current_Adaper< Ty_HH_FT_GH_CUR<TySineCurrent> > Ty_HH_FT_GH_sine;


struct Ty_HH_GH_cont_syn
  :public Ty_Neuron_Dym_Base
{
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
  static const int id_h     = 1;
  static const int id_m     = 2;
  static const int id_n     = 3;
  static const int id_gE    = 4;
  static const int id_gI    = 5;
  static const int id_gE_s1 = 6;
  static const int id_gI_s1 = 7;
  static const int id_gEInject = id_gE_s1;
  static const int id_gIInject = id_gI_s1;

  double Get_V_threshold() const override {return V_threshold;};
  int Get_id_gEInject() const override {return id_gEInject;}
  int Get_id_gIInject() const override {return id_gIInject;}
  int Get_n_dym_vars() const override {return n_var;}
  int Get_id_V() const override {return id_V;}
  int Get_id_gE() const override {return id_gE;}
  int Get_id_gI() const override {return id_gI;}
  void NextStepSingleNeuronQuiet(
    double *dym_val,
    double &t_in_refractory,
    double &spike_time_local,
    double dt_local) const override
  {}
  void VoltHandReset(double *dym_val) const override
  {}
};

#endif
