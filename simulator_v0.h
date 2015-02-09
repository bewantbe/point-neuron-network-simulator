
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

