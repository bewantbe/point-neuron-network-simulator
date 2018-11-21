// Read input options and run the simulation

/*
TODO:
   * correctness check.
   * add an interface to convert the parameters from mV(PSP) to strength. (so easier to read)
   * divide NDEBUG in to multiple levels.
   * speed/accuracy benchmark.
   * HH model: use the time of spike top instead of threshold for synaptic interaction.
   * speed up GetDv() in HH model. vectorize the exp() for further speed up.
   * use sorted spike detection and partial update to speed up?
   * merge t_in_ref to dym_val?
   * possibility to add synaptic delay? Partly done.
   * add to Poisson event queue?
   * isomerisom of neurons in a network?
*/



#ifndef DEBUG
#define NDEBUG  // disable assert() and disable checks in Eigen
#endif

#include "common_header.h"
#include "legacy_lib.h"

#include "single_neuron_dynamics.h"
#include "neuron_population.h"
#include "simulator_exact_order.h"

#include "neuron_population_cont_synaptic.h"
#include "simulator_cont_synaptic.h"
//FILE *fdbg;
#include "simulator_if_jump.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef std::vector<unsigned int> VecUInt;

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

size_t quiet_step_call_counter = 0;

#include <stdarg.h>

int tmp_dbg_printf(const char *format, ...)
{
  static long cnt = 0;
  cnt++;
  if (cnt < (long)(10000l * 25 * 78*32)) return 0;
  va_list ap;
  va_start(ap, format);
  int ret = vfprintf(stdout, format, ap);
  va_end(ap);
  return ret;
}

// Define wall clock
#ifdef _WINDOWS_
#include <time.h>
clock_t tic()
{ return clock(); }
double toc(const clock_t &t_begin)
{ return (clock()-t_begin)/(1.0*CLOCKS_PER_SEC); }
#else
#include <sys/time.h>
timeval tic()
{
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv;
}
double toc(const timeval &t0)
{
  timeval t1;
  gettimeofday(&t1, NULL);
  return (t1.tv_sec - t0.tv_sec)+1e-6*(t1.tv_usec - t0.tv_usec);
}
#endif

void ForceSpikeOnList(NeuronPopulationBase *p_neu_pop,
                      const TySpikeEventVec &force_spike_list, double t, double dt)
{
  static size_t id_next_spike = 0;
  if (p_neu_pop == nullptr) {
    id_next_spike = 0;
    return;
  }
  if (id_next_spike >=force_spike_list.size())
    return;

  const double t_fire = force_spike_list[id_next_spike].time;
  if (t <= t_fire && t_fire < t+dt) {
    p_neu_pop->ForceReset(force_spike_list[id_next_spike].id);
    id_next_spike++;
  }
}

int MainLoop(const po::variables_map &vm)
{
  auto t_begin = tic();
  if (!vm.count("neuron-model")) {
    cerr << "Error: Use \"--neuron-model\" to specify a neuron model. See --help.\n";
    return -1;
  }
  const std::string &str_nm = vm["neuron-model"].as<std::string>();

  bool b_verbose = false;
  if (vm.count("verbose")) {
    b_verbose = true;
  }
  
  bool b_verbose_echo = false;
  if (vm.count("verbose-echo")) {
    b_verbose_echo = true;
  }

  // Neuron models.
  enum EnumNeuronModel {
    LIF_G, LIF_GH, HH_G, HH_GH, HH_PT_GH, HH_FT_GH,
    HH_G_sine, HH_GH_sine, HH_PT_GH_sine, HH_FT_GH_sine,
    HH_G_extI, HH_GH_extI,
    HH_GH_cont_syn, IF_jump, Hawkes_GH,
    DIF_single_GH, DIF_GH
  };
  EnumNeuronModel enum_neuron_model;
  if (str_nm ==            "LIF-G") {
    enum_neuron_model =     LIF_G;
  } else if (str_nm ==     "LIF-GH") {
    enum_neuron_model =     LIF_GH;
  } else if (str_nm ==     "HH-G") {
    enum_neuron_model =     HH_G;
  } else if (str_nm ==     "HH-GH") {
    enum_neuron_model =     HH_GH;
  } else if (str_nm ==     "HH-PT-GH") {
    enum_neuron_model =     HH_PT_GH;
  } else if (str_nm ==     "HH-FT-GH") {
    enum_neuron_model =     HH_FT_GH;
  } else if (str_nm ==     "HH-G-sine") {
    enum_neuron_model =     HH_G_sine;
  } else if (str_nm ==     "HH-GH-sine") {
    enum_neuron_model =     HH_GH_sine;
  } else if (str_nm ==     "HH-PT-GH-sine") {
    enum_neuron_model =     HH_PT_GH_sine;
  } else if (str_nm ==     "HH-FT-GH-sine") {
    enum_neuron_model =     HH_FT_GH_sine;
  } else if (str_nm ==     "HH-G-extI") {
    enum_neuron_model =     HH_G_extI;
  } else if (str_nm ==     "HH-GH-extI") {
    enum_neuron_model =     HH_GH_extI;
  } else if (str_nm ==     "HH-GH-cont-syn") {
    enum_neuron_model =     HH_GH_cont_syn;
  } else if (str_nm ==     "IF-jump") {
    enum_neuron_model =     IF_jump;
  } else if (str_nm ==     "Hawkes-GH") {
    enum_neuron_model =     Hawkes_GH;
  } else if (str_nm ==     "DIF-single-GH") {
	enum_neuron_model =		DIF_single_GH;
  } else if (str_nm ==     "DIF-GH"){
	  enum_neuron_model = DIF_GH;
  } else {
    cerr << "Unrecognized neuron model. See --help.\n";
    return -1;
  }

  // Set neural parameters
  TyNeuronalParams pm(vm["nE"].as<unsigned int>(), vm["nI"].as<unsigned int>());

  // Set poisson input rate and strength
  auto f_parse_array = [&vm](const char *opt_name, std::vector<double> &to_arr)
  {
    if (vm.count(opt_name)) {
      const auto &v = vm[opt_name].as< std::vector<double> >();
      if (v.size() == 1) {
        std::fill(to_arr.begin(), to_arr.end(), v[0]);
      } else {
        size_t upbd = std::min(v.size(), to_arr.size());
        for (size_t i = 0; i < upbd; i++) {
          to_arr[i] = v[i];
        }
      }
    }
    //printf("Arr %s\n", opt_name);
    //for (size_t i = 0; i < to_arr.size(); i++) {
      //printf("  arr[%lu] = %.17g\n", i, to_arr[i]);
    //}
  };

  f_parse_array("pr", pm.arr_pr);
  f_parse_array("ps", pm.arr_ps);
  f_parse_array("pri", pm.arr_pri);
  f_parse_array("psi", pm.arr_psi);

  pm.scee = vm["scee"].as<double>();
  pm.scie = vm["scie"].as<double>();
  pm.scei = vm["scei"].as<double>();
  pm.scii = vm["scii"].as<double>();

  // Initialize net before initialize delay-net (in neuron population)
  if (vm.count("net")) {
    std::string name_net = vm["net"].as<std::string>();
    bool is_sparse = vm.count("sparse-net") > 0;
    FillNetFromPath(pm, name_net, is_sparse);
  } else {
    cout << "You must specify a network. (--net)" << endl;
    return -1;
  }

  // Only for DIF
  // Parameters about --synaptic-alpha-path
  // Initialize alpha_coef(ficient) before DIF-related models contruction
  pm.DIF_flag = false;
  if (enum_neuron_model == DIF_GH) {
	  pm.DIF_flag = true;
	  if (vm.count("synaptic-alpha-path")) {
		  std::string name_coef = vm["synaptic-alpha-path"].as<std::string>();
		  bool is_sparse = vm.count("sparse-net") > 0;
		  InitAlphaCoeffFromPath(pm, name_coef, is_sparse);
	  }
	  else {
		  cout << "You must specify synaptic alpha(s). (--synaptic-alpha-path)" << endl;
		  return -1;
	  }
  }

  // Set neuron population type
  NeuronPopulationBase * p_neu_pop = nullptr;
  if (vm.count("synaptic-net-delay")) {
    switch (enum_neuron_model) {
      case LIF_G:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_LIF_G>(pm);
        break;
      case LIF_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_LIF_GH>(pm);
        break;
      case HH_G:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_HH_G>(pm);
        break;
      case HH_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_HH_GH>(pm);
        break;
      case HH_PT_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_HH_PT_GH>(pm);
        break;
      case HH_FT_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractNetDelay<Ty_HH_FT_GH>(pm);
        break;
      case IF_jump:
        p_neu_pop = new IFJumpPopulation(pm);
        break;
      default:
        cerr << "The neuron type with delay net is not supported yet.\n";
        return -1;
    };
    std::string path = vm["synaptic-net-delay"].as<std::string>();
    try {
      SparseMat dn = ReadNetDelay(path, pm.net);
      p_neu_pop->SetSynapticDelayNet(dn);
    } catch (const char*) {
      if (vm.count("synaptic-delay")) {
        // if ReadNetDelay() somehow failed, we can still do this.
        cout << "Warning: Using constant delay.\n";
        p_neu_pop->SetSynapticDelay(vm["synaptic-delay"].as<double>());
      }
    }
  } else if (vm.count("synaptic-delay")) {
    switch (enum_neuron_model) {
      case LIF_G:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_LIF_G>(pm);
        break;
      case LIF_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_LIF_GH>(pm);
        break;
      case HH_G:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_HH_G>(pm);
        break;
      case HH_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_HH_GH>(pm);
        break;
      case HH_PT_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_HH_PT_GH>(pm);
        break;
      case HH_FT_GH:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelay<Ty_HH_FT_GH>(pm);
        break;
      case HH_G_sine:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelaySine<Ty_HH_G_sine>(pm);
        break;
      case HH_GH_sine:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelaySine<Ty_HH_GH_sine>(pm);
        break;
      case HH_PT_GH_sine:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelaySine<Ty_HH_PT_GH_sine>(pm);
        break;
      case HH_FT_GH_sine:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelaySine<Ty_HH_FT_GH_sine>(pm);
        break;
      case HH_G_extI:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelayExtI<Ty_HH_G_extI>(pm);
        break;
      case HH_GH_extI:
        p_neu_pop = new
          NeuronPopulationDeltaInteractConstantDelayExtI<Ty_HH_GH_extI>(pm);
        break;
      case HH_GH_cont_syn:
        cerr << "Delay for HH_GH_cont_syn is not supported yet.\n";
        return -1;
      case IF_jump:
        p_neu_pop = new IFJumpPopulation(pm);
        return -1;
      case Hawkes_GH:
        cerr << "Delay for Hawkes-GH is not supported yet.\n";
        return -1;
      default:
        throw "not supported combination of model and network.";
    }
    p_neu_pop->SetSynapticDelay(vm["synaptic-delay"].as<double>());
  } else if (vm.count("tau-g-path") || vm.count("neuron-const-path")) {
    switch (enum_neuron_model) {
      case LIF_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractHeterogeneous<Ty_LIF_GH>(pm);
        break;
      case DIF_single_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractHeterogeneous<Ty_DIF_single_GH>(pm);
        break;
      case HH_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractHeterogeneous<Ty_HH_GH>(pm);
        break;
      default:
        throw "tau-g-path: only LIF-GH, DIF_single_GH and HH-GH are supported";
    }
    if (vm.count("tau-g-path")) {
      TyNeuDymParam dym_param;
      FillTauG(dym_param, vm["tau-g-path"].as<std::string>().c_str());
      p_neu_pop->SetRisingFallingTau(dym_param);
    }
    if (vm.count("neuron-const-path")) {
      TyArrVals v;
      FillVecFromFile(v, vm["neuron-const-path"].as<std::string>().c_str());
      p_neu_pop->SetAllDymParam(v);
    }
  } else {
    switch (enum_neuron_model) {
      case LIF_G:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_LIF_G>(pm);
        break;
      case LIF_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_LIF_GH>(pm);
        break;
      case HH_G:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_HH_G>(pm);
        break;
      case HH_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_HH_GH>(pm);
        break;
      case HH_PT_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_HH_PT_GH>(pm);
        break;
      case HH_FT_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_HH_FT_GH>(pm);
        break;
      case HH_G_sine:
        p_neu_pop = new NeuronPopulationDeltaInteractSine<Ty_HH_G_sine>(pm);
        break;
      case HH_GH_sine:
        p_neu_pop = new NeuronPopulationDeltaInteractSine<Ty_HH_GH_sine>(pm);
        break;
      case HH_PT_GH_sine:
        p_neu_pop = new NeuronPopulationDeltaInteractSine<Ty_HH_PT_GH_sine>(pm);
        break;
      case HH_FT_GH_sine:
        p_neu_pop = new NeuronPopulationDeltaInteractSine<Ty_HH_FT_GH_sine>(pm);
        break;
      case HH_G_extI:
        p_neu_pop = new NeuronPopulationDeltaInteractExtI<Ty_HH_G_extI>(pm);
        break;
      case HH_GH_extI:
        p_neu_pop = new NeuronPopulationDeltaInteractExtI<Ty_HH_GH_extI>(pm);
        break;
      case HH_GH_cont_syn:
        p_neu_pop = new NeuronPopulationContSyn(pm);
        break;
      case IF_jump:
        p_neu_pop = new IFJumpPopulation(pm);
        break;
      case Hawkes_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_Hawkes_GH>(pm);
        break;
      case DIF_single_GH:
        p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_DIF_single_GH>(pm);
        break;
	  case DIF_GH:
		p_neu_pop = new NeuronPopulationDendriticDeltaInteractTemplate<Ty_DIF_GH>(pm);
    }
  }
  // Get basic single neuron info
  const Ty_Neuron_Dym_Base *p_neuron_model = p_neu_pop->GetNeuronModel();

  double e_t   = vm["t"].as<double>();
  double e_dt  = vm["dt"].as<double>();
  double e_stv = e_dt;
  if (vm.count("stv")) {
    e_stv = vm["stv"].as<double>();
  }
  if (enum_neuron_model == IF_jump) {
    // There is no point to use finer dt in Poisson driven IF-jump model.
    if (e_stv > e_dt && e_dt < 1) {
      cout << "Warning: dt (" << e_dt << " ms) < stv (" << e_stv << " ms)";
      cout << " in IF-jump model, which is pointless. Use a larger dt please.\n";
    }
  }
  if (!(e_t > 0) || !(e_dt > 0) || !(e_stv > 0)) {
    cerr << "Only support positive time!" << endl;
    return 2;
  }
  if (!(e_dt<=e_stv && e_stv<=e_t)) {
    cerr << "Must dt <= stv <= t !" << endl;
    return 2;
  }
  if (e_stv / floor(e_stv/e_dt + 0.5) != e_dt) {  // stv is not a multiple of dt
    cerr << "stv must be a multiple of dt !" << endl;
    return 2;
  }
  if (e_stv / e_dt > INT_MAX) {
    cerr << "stv / dt > INT_MAX ! (dt too small or stv too large)" << endl;
    return 2;
  }

  // Parameters about sine current input
  if (enum_neuron_model == HH_GH_sine
      || enum_neuron_model == HH_G_sine
      || enum_neuron_model ==  HH_FT_GH_sine) {
    auto p_neu_pop_sine = dynamic_cast<
        NeuronPopulationBaseSine *>(p_neu_pop);
    if (vm.count("current-sine-amp")) {
      p_neu_pop_sine->SetSineAmplitude(
          vm["current-sine-amp"].as<double>());
    }
    if (vm.count("current-sine-freq")) {
      p_neu_pop_sine->SetSineFrequency(
          vm["current-sine-freq"].as<double>());
    }
  }

  // Parameters about -extI current input
  if (enum_neuron_model == HH_G_extI
      ||enum_neuron_model == HH_GH_extI) {
    auto p_neu_pop_extI = dynamic_cast<
      NeuronPopulationBaseExtI *>(p_neu_pop);
    if (vm.count("extI")) {
      p_neu_pop_extI->SetExtI(
          vm["extI"].as<double>());
    }
  }

  // Set simulator for the neuronal network.
  NeuronSimulatorBase *p_neu_simu = nullptr;
  std::string str_simu_method = "";
  if (vm.count("simulation-method")) {
    str_simu_method = vm["simulation-method"].as<std::string>();
  }
  double t_warming_up = vm["t-warming-up"].as<double>();
  // Align t_warming_up to dt, so that t_warming_up is easier reachable and
  // also time points for volt are multiples of dt.
  t_warming_up = floor(t_warming_up / e_dt) * e_dt;
  double t0 = 0; // TODO: t0 may used to account for non-aligned t_warming_up
  if (str_simu_method.size() > 0 && str_simu_method != "auto") {
    str_simu_method = vm["simulation-method"].as<std::string>();
    if (str_simu_method == "simple") {
      p_neu_simu = new NeuronSimulatorSimple(pm, e_dt, t0);
    } else if (str_simu_method == "SSC") {  // Spike-Spike-Correction
      p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt, t0);
    } else if (str_simu_method == "SSC-Sparse") {
      p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse(pm, e_dt, t0);
    } else if (str_simu_method == "SSC-Sparse2") {
      p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse2(pm, e_dt, t0);
    } else if (str_simu_method == "big-net-delay") {
      p_neu_simu = new NeuronSimulatorBigNetDelay(pm, e_dt, t0);
    } else if (str_simu_method == "big-delay") {
      p_neu_simu = new NeuronSimulatorBigDelay(pm, e_dt, t0);
    } else if (str_simu_method == "cont-syn") {
      p_neu_simu = new NeuronSimulatorCont(pm, e_dt, t0);
    } else if (str_simu_method == "IF-jump") {
      p_neu_simu = new IFJumpSimulator(pm, e_dt, t0);
    } else if (str_simu_method == "IF-jump-delay") {
      p_neu_simu = new IFJumpDelayNetSimulator(pm, e_dt, t0);
    } else {
      cerr << "No this simulation method:\"" << str_simu_method << "\"\n";
      return -2;
    }
  } else {
    // Default simulator
    if (HH_GH_cont_syn == enum_neuron_model) {
      str_simu_method = "cont-syn";
      p_neu_simu = new NeuronSimulatorCont(pm, e_dt, t0);
    } else if (IF_jump == enum_neuron_model) {
      if (vm.count("synaptic-net-delay") || vm.count("synaptic-delay")) {
        str_simu_method = "IF-jump-delay";
        p_neu_simu = new IFJumpDelayNetSimulator(pm, e_dt, t0);
      } else {
        str_simu_method = "IF-jump";
        p_neu_simu = new IFJumpSimulator(pm, e_dt, t0);
      }
    } else if (vm.count("synaptic-net-delay")) {
      str_simu_method = "big-net-delay";
      p_neu_simu = new NeuronSimulatorBigNetDelay(pm, e_dt, t0);
    } else if (vm.count("synaptic-delay")) {
      str_simu_method = "big-delay";
      p_neu_simu = new NeuronSimulatorBigDelay(pm, e_dt, t0);
    } else {
      str_simu_method = "SSC";
      p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt, t0);
    }
    //cout << "Simulator: " << str_simu_method << "\n";  // Debug info
  }
  // Check inconsistant pair
  if (vm.count("synaptic-delay") && 
      (str_simu_method != "big-delay" && str_simu_method != "big-net-delay")) {
    cerr << "You select a delayed synaptic network, but use a non-compatible simulator. Try --simulation-method big-delay\n";
    return -1;
  }
  if ((enum_neuron_model == HH_GH_cont_syn)
      ^ (str_simu_method == "cont-syn")) {
    cerr << "Error: Neuron Model \"HH_GH_cont_syn\" must pair with simulator \"cont-syn\".\n";
    return -1;
  }
  if ((enum_neuron_model == IF_jump)
      ^ (str_simu_method == "IF-jump" || str_simu_method == "IF-jump-delay")) {
    cerr << "Error: Neuron Model \"IF-jump\" must pair with simulator \"IF-jump\" or \"IF-jump-delay\".\n";
    return -1;
  }
  if (enum_neuron_model == Hawkes_GH && str_simu_method != "simple") {
    cerr << "Currently Hawkes-GH model only support \"simple\" simulator.\n";
    return -1;
  }
  
  if (vm.count("set-threshold")) {
    p_neu_pop->SetThreshold(vm["set-threshold"].as<double>());
  }
  if (vm.count("no-threshold")) {
    p_neu_pop->DisableThreshold();
  }

  bool b_init_loaded = false;
  if (vm.count("initial-state-path")) {
    int rt = FillNeuStateFromFile(p_neu_pop->GetDymState(),
                 vm["initial-state-path"].as<std::string>().c_str());
    if (rt == 0) {
      b_init_loaded = true;
      if (b_verbose_echo) cout << "Initial state loaded." << endl;
    }
    //p_neu_simu->SaneTestVolt();
  }
  if (!b_init_loaded) {
    // Fill with default values
    auto &ds = p_neu_pop->GetDymState().dym_vals;
    int n_var = p_neuron_model->Get_n_dym_vars();
    for (int i = 0; i < p_neu_pop->n_neurons(); i++) {
      const double *v = p_neuron_model->Get_dym_default_val();
      memcpy(ds.data()+i*n_var, v, n_var*sizeof(double));
    }
    b_init_loaded = true;
  }
//    const auto &ds = p_neu_pop->GetDymState().dym_vals;
//    for (int i = 0; i < ds.rows(); i++) {
//      printf("%.16e ", ds(i, 0));
//      printf("%.16e ", ds(i, 1));
//      printf("%.16e ", ds(i, 2));
//      printf("%.16e ", ds(i, 3));
//      printf("%.16e ", ds(i, 4));
//      printf("%.16e ", ds(i, 5));
//      printf("%.16e ", ds(i, 6));
//      printf("%.16e ", ds(i, 7));
//      printf("\n");
//    }

  if (vm.count("input-event-path")) {
    FillPoissonEventsFromFile(p_neu_simu->Get_poisson_time_vec(),
                              vm["input-event-path"].as<std::string>().c_str(),
                              p_neu_pop->GetNeuronalParamsPtr()->arr_ps);
  }

  auto fout_try_open = [&vm](const char * const st_id, std::ofstream &fs)
    -> bool {
      if (vm.count(st_id)) {
        std::string st_path = vm[st_id].as<std::string>();
        CheckDirAndCreate( st_path );
        fs.open( st_path, std::ios::binary );
        if (!fs) {
          cerr << "Error: Failed to open file \"" << st_path << "\" for output." << "\n";
          throw "Failed to open file";
        }
        return true;
      } else {
        return false;
      }
    };

  // Try open the files for data dump
  std::ofstream fout_volt;
  bool b_output_volt = fout_try_open("volt-path", fout_volt);

  std::ofstream fout_G;    // G: conductance
  bool b_output_G = fout_try_open("conductance-path", fout_G);

  bool b_output_extra_G_DIF = (enum_neuron_model == DIF_GH);

  std::ofstream fout_HH_gate;    // File for saving gating variables: h, m, n
  bool b_output_HH_gate = false;
  if (str_nm.find("HH") != std::string::npos) {
    b_output_HH_gate = fout_try_open("ion-gate-path", fout_HH_gate);
  }

  std::ofstream fout_ras;
  bool output_ras = fout_try_open("ras-path", fout_ras);
  if (output_ras) {
    fout_ras.precision(17);
  }
  
  std::ofstream fout_poisson;
  bool b_output_poisson_input = fout_try_open("output-poisson-path", fout_poisson);
  if (b_output_poisson_input) {
    fout_poisson.precision(17);
  }

  // 
  TySpikeEventVec force_spike_list;
  if (vm.count("force-spike-list")) {
    ReadSpikeList(force_spike_list, 
      vm["force-spike-list"].as<std::string>().c_str());
  }

  // Set Time_refractory
  if (vm.count("refractory-time")) {
    p_neu_pop->SetRefractoryTime(vm["refractory-time"].as<double>());
  }

  std::vector<size_t> vec_n_spike(p_neu_pop->n_neurons());  // count the number of spikes
  TySpikeEventVec ras;                            // record spike raster
  int n_dt_in_stv = int(e_stv / e_dt + 0.1);
  int count_n_dt_in_stv = n_dt_in_stv;
  size_t n_step = (size_t)((e_t + t_warming_up) / e_dt);
  
  // Show input parameters.
  if (b_verbose_echo) {
    printf("Model: \"%s\"\n", str_nm.c_str());
    const TyNeuronalParams &cpm = *p_neu_pop->GetNeuronalParamsPtr();
    printf("Number of neurons: %d E + %d I\n", cpm.n_E, cpm.n_I);
    int show_n = std::min(10, cpm.n_total());
    unsigned long n_nz = cpm.net.nonZeros();
    unsigned long n_nm = (unsigned long)(cpm.net.rows()-1) * cpm.net.cols();
    printf("Network: number of edges = %lu/%lu (%.2g %%)\n", n_nz, n_nm,
           100.0*n_nz/n_nm);
    for (int i = 0; i < show_n; i++) {
      for (int j = 0; j < show_n; j++) {
        printf("% 8.2g", cpm.net.coeff(i, j));
      }
      printf("\n");
    }
    printf("pr :\n");
    for (int i = 0; i < show_n; i++) {
      printf("% 8.2g", cpm.arr_pr[i]);
    }
    printf("\n");
    printf("ps :\n");
    for (int i = 0; i < show_n; i++) {
      printf("% 8.2g", cpm.arr_ps[i]);
    }
    printf("\n");
    auto any = [](std::vector<double> v) {
      return std::any_of(v.cbegin(), v.cend(),
        [](double i){ return i != 0; }
      );
    };
    if (any(cpm.arr_psi)) {
      printf("pri:\n");
      for (int i = 0; i < show_n; i++) {
        printf("% 8.2g", cpm.arr_pri[i]);
      }
      printf("\n");
      printf("psi:\n");
      for (int i = 0; i < show_n; i++) {
        printf("% 8.2g", cpm.arr_psi[i]);
      }
      printf("\n");
    }
    printf("\n     scee     scie     scei     scii\n");
    printf(" %8.3g %8.3g %8.3g %8.3g\n",
           cpm.scee, cpm.scie, cpm.scei, cpm.scii);
    printf("\nSimulator: \"%s\"\n", str_simu_method.c_str());
    if (t_warming_up > 0)
      printf("  t = %.2g (warming-up) + %.2g ms, dt = %.2g ms, stv = %.2g ms (%d dt)\n",
             t_warming_up, e_t, e_dt, e_stv, n_dt_in_stv);
    else
      printf("  t = %.2g ms, dt = %.2g ms, stv = %.2g ms (%d dt)\n",
             e_t, e_dt, e_stv, n_dt_in_stv);
    if (vm.count("input-event-path")) {
      TyPoissonTimeVec &pv = p_neu_simu->Get_poisson_time_vec();
      unsigned long n_events = 0;
      for (size_t j = 0; j < pv.size(); j++) {
        n_events += pv[j].size() - 1;
      }
      printf("  Imported events: %lu\n", n_events);
    }
    if (b_verbose) printf("\n");
  }

  // Show progress.
  if (b_verbose) {
    printf("  Initialization     : %3.3f s\n", toc(t_begin));
    fflush(stdout);
    printf("  Simulation         : ");
  }


  /*                             TODOYWS: check
	 example of the case storing conductance of DIF-GH model (n neuron, [t/dt] = m):
	 {t = 1*dt, 1st neu} {t = 1*dt, 2ed neu} ... {t = 1*dt, nst neu} 
	 {t = 2*dt, 1st neu} {t = 2*dt, 2ed neu} ... {t = 2*dt, nst neu}
	 ......
	 {t = m*dt, 1st neu} {t = m*dt, 2ed neu} ... {t = m*dt, nst neu}
	 where
	 {t = i*dt, jth ner} means: data of jth neu at time i*dt,
	 'data' includes gEP gIP gEPs gIPs  gE(1) gE(2) ... gE(n)  gI(1) gI(2) ... gI(n), refer to @@Ty_DIF_GH_core for more information
   */
  // Function for save data to file
  auto func_save_dym_values = [
    b_output_volt, &fout_volt,
    b_output_G, b_output_extra_G_DIF, &fout_G,     // EIDTED YWS, to enable output extra conductance data for DIF model
    b_output_HH_gate, &fout_HH_gate,
    p_neuron_model]
    (const NeuronPopulationBase &neu_pop) {
    if (b_output_volt) {
      std::vector<double> v(neu_pop.n_neurons());
      int id_V = p_neuron_model->Get_id_V();
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        v[j] = neu_pop.GetDymState().dym_vals(j, id_V);
      }
      fout_volt.write((char*)v.data(), neu_pop.n_neurons()*sizeof(double));
    }
    if (b_output_G && ! b_output_extra_G_DIF) {
      int id_gE = p_neuron_model->Get_id_gE();
      int id_gI = p_neuron_model->Get_id_gI();
      std::vector<double> v(2*neu_pop.n_neurons());
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        v[2*j  ] = neu_pop.GetDymState().dym_vals(j, id_gE);
        v[2*j+1] = neu_pop.GetDymState().dym_vals(j, id_gI);
      }
      fout_G.write((char*)v.data(), 2*neu_pop.n_neurons()*sizeof(double));
    }
	if (b_output_G && b_output_extra_G_DIF) {							// output conductence for DIF-GH model

		//int n_var = (2 + 2 * neu_pop.n_neurons()) * neu_pop.n_neurons();
		//std::vector<double> v(n_var*neu_pop.n_neurons());
		//for (int j = 0; j < neu_pop.n_neurons(); j++) {
			//v[n_var*j + 0] = neu_pop.GetDymState().dym_vals(j, 1);      // id_gEPoisson == 1
			//v[n_var*j + 1] = neu_pop.GetDymState().dym_vals(j, 2);		// id_gIPoisson == 2
			//for (int i = 2; i < n_var; i++) {
			//	v[n_var*j + i] = neu_pop.GetDymState().dym_vals(j, i+3); // we skipped '2' varibles, namely gEPs, gIPs
			//}

		//}

		
		   // the following part of code output all dynamic varibles, including gEP gEI gE(i) gI(i) gEs(i) gIs(i), excluding V
		   // useful when debug
		int n_var = (p_neuron_model->Get_n_dym_vars() - 1 + 4 * neu_pop.n_neurons()); // FIXED YWS 2018/11/15  * neu_pop.n_neurons();
		std::vector<double> v(n_var*neu_pop.n_neurons());
		for (int j = 0; j < neu_pop.n_neurons(); j++) {					// I think we can just write the whole dym_val[]
			for (int i = 0; i < n_var; i++) {
				v[n_var*j + i] = neu_pop.GetDymState().dym_vals(j, i+1);
			}
		}
		
		fout_G.write((char*)v.data(), n_var*neu_pop.n_neurons()*sizeof(double));

	}
    if (b_output_HH_gate) {
      int id_gating = p_neuron_model->Get_id_V() + 1;
      int n_gating = ((Ty_HH_GH*)p_neuron_model)->n_var_soma - 1;
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        const double *p = & neu_pop.GetDymState().dym_vals(j, id_gating);
        fout_HH_gate.write((char*)p, sizeof(double)*n_gating);
      }
    }
  };

  if (vm.count("output-first-data-point")) {
    func_save_dym_values(*p_neu_pop);
  }

  // Main loop
  t_begin = tic();
  int progress_percent = 0;
  int str_len = 0;
  for (size_t i = 0; i < n_step; i++) {
    ForceSpikeOnList(p_neu_pop, force_spike_list, p_neu_simu->GetT(), e_dt);
    if (b_output_poisson_input) {
      SavePoissonInput(fout_poisson, p_neu_simu->Get_poisson_time_vec(),
          p_neu_simu->GetT()+e_dt);
    }
    // Advance one dt step
    p_neu_simu->NextDt(p_neu_pop, ras, vec_n_spike);
    
    // Show progress
    if (b_verbose && (i+1.0)/n_step >= (progress_percent+1.0)/100) {
      progress_percent = (int)((i+1.0)/n_step * 100);
      // For matlab compatible (otherwise just printf('\r') will be fine)
      for (int j=0; j<str_len; j++) printf("\b");
      str_len = printf("%3.3f s (%d%%)",
                       toc(t_begin), progress_percent);
      fflush(stdout);
    }

    // Discard warming-up data
    // Use (i+1)*e_dt instead of p_neu_simu->GetT() to make "==" possible.
    if ((i+1)*e_dt <= t_warming_up) {
      ras.clear();
      std::fill(vec_n_spike.begin(), vec_n_spike.end(), 0);
      continue;
    }

    if (output_ras) {
      for (size_t j = 0; j < ras.size(); j++) {
        fout_ras << ras[j].id + 1 << '\t' << ras[j].time - t_warming_up << '\n';
      }
    }
    ras.clear();
    // output dynamical variable(s) every n_dt_in_stv
    if (--count_n_dt_in_stv > 0)
      continue;
    else
      count_n_dt_in_stv = n_dt_in_stv;
    func_save_dym_values(*p_neu_pop);
  }
  if (b_verbose) {
    printf("\n");
    fflush(stdout);
  }
  
  if (b_verbose_echo) {
    if (b_verbose) printf("\n");
    const TyNeuronalParams &cpm = *p_neu_pop->GetNeuronalParamsPtr();
    int show_n = std::min(10, cpm.n_total());
    printf("Number of spikes:\n");
    for (int i = 0; i < show_n; i++) {
      printf(" %7lu", vec_n_spike[i]);
    }
    printf("\n");
    
    printf("# of calls to NoInteractDt(): %lu\n", quiet_step_call_counter);
  }

  if (vm.count("isi-path")) {
    std::ofstream fout_isi( vm["isi-path"].as<std::string>() );
    fout_isi.precision(17);
    for (int j = 0; j < p_neu_pop->n_neurons(); j++) {
      fout_isi << e_t / vec_n_spike[j] << '\t';
    }
  }

//  delete p_neu_simu;

  return 0;
}

int main(int argc, char *argv[])
{
//  fdbg = fopen("dbg.txt", "w");
  // Declare the supported options.
  po::options_description desc("Options");
  // http://stackoverflow.com/questions/3621181/short-options-only-in-boostprogram-options
  desc.add_options()
      ("neuron-model",  po::value<std::string>(),
       "One of LIF-G, LIF-GH, HH-G, HH-GH, HH-PT-GH, HH-FT-GH, HH-G-sine, HH-GH-sine, HH-PT-GH-sine, HH-FT-GH-sine, HH-G-extI, HH-GH-extI, HH-GH-cont-syn, IF-jump, Hawkes-GH, DIF-single-GH, DIF-GH.")
      ("simulation-method",  po::value<std::string>(),
       "One of simple, SSC, SSC-Sparse, SSC-Sparse2, big-delay, big-net-delay, cont-syn, IF-jump, IF-jump-delay, auto. Some combinations of neuron model and simulator are mutually exclusive, hence not allowed. If not specify, a suitable simulator will be choosen automatically.")
      ("help,h",
       "Produce help message.")
      ("verbose,v",
       "Show progress.")
      ("verbose-echo",
       "Show input parameters.")
      ("t",    po::value<double>()->default_value(1e3),
       "Simulation time, in ms.")
      ("dt",   po::value<double>()->default_value(1.0/32),
       "Simulation delta t (dt, time step), in ms.")
      ("stv",  po::value<double>(),
       "Output sampling interval. Must be multiples of dt. Default set to dt.")
      ("t-warming-up", po::value<double>()->default_value(0.0),
       "The initial time period to be discarded, in ms.")
      ("nE",   po::value<unsigned int>()->default_value(1),
       "Number of excitatory neurons.")
      ("nI",   po::value<unsigned int>()->default_value(0),
       "Number of inhibitory neurons.")
      ("net",  po::value<std::string>(),
       "Network file path. Use - for full net.")
      ("sparse-net",
       "The network file is in sparse format, i.e. i j value.")
      ("scee", po::value<double>()->default_value(0.0),
       "Cortical strength E->E (unscaled).")
      ("scie", po::value<double>()->default_value(0.0),
       "Cortical strength E->I (unscaled).")
      ("scei", po::value<double>()->default_value(0.0),
       "Cortical strength I->E (unscaled).")
      ("scii", po::value<double>()->default_value(0.0),
       "Cortical strength I->I (unscaled).")
      ("pr",   po::value< std::vector<double> >()->multitoken(),
       "Poisson input rate for Excitatory type, in 1/ms.")
      ("ps",   po::value< std::vector<double> >()->multitoken(),
       "Poisson input strength for Excitatory type (unscaled).")
      ("pri",  po::value< std::vector<double> >()->multitoken(),
       "Poisson input rate for inhibitory type, in 1/ms.")
      ("psi",  po::value< std::vector<double> >()->multitoken(),
       "Poisson input strength, inhibitory type (unscaled).")
      ("seed", po::value< std::vector<std::string> >()->multitoken(),
       "Random seed for Poisson events. One or several unsigned integers (0 ~ 2^32-1), or 'auto'.")
      ("set-threshold", po::value<double>(),
       "Set the threshold to non-standard value.")
      ("no-threshold",
       "Disable the threshold detection and voltage reset. Note that HH based neurons will still fire.")
      ("current-sine-amp", po::value<double>(),
       "Set the current input sine amplitude.")
      ("current-sine-freq", po::value<double>(),
       "Set the current input sine frequency, in kHz.")
      ("extI", po::value<double>(),
       "Set the custom current input extI.")
      ("synaptic-delay", po::value<double>(),
       "Set a synaptic delay for the network.")
      ("synaptic-net-delay", po::value<std::string>(),
       "Set the synaptic delays for the network. Path to the matrix text file.")
      ("refractory-time", po::value<double>(),
       "Set refractory time for LIF type models.")
      ("volt-path,o",      po::value<std::string>(),
       "Output volt to path. In raw binary format.")
      ("ras-path",         po::value<std::string>(),
       "Output spike events to path.")
      ("isi-path",         po::value<std::string>(),
       "Output mean InterSpike Interval to path.")
      ("conductance-path", po::value<std::string>(),
       "Output conductance to path")
      ("ion-gate-path", po::value<std::string>(),
       "Output gating variables to path. In raw binary format.")
      ("output-first-data-point",
       "Also output initial condition in volt-path and conductance-path.")
      ("output-poisson-path", po::value<std::string>(),
       "Output generated poisson events to the path")
      ("tau-g-path", po::value<std::string>(),
       "Set rising and falling time constants for E and I conductance for each neuron from path.")
      ("neuron-const-path", po::value<std::string>(),
       "Set dynamical constants for each neuron from path.")
      ("initial-state-path", po::value<std::string>(),
       "Read initial state from path.")
      ("input-event-path", po::value<std::string>(),
       "Read input event from path.")
      ("parameter-path", po::value<std::string>(),
       "Read parameters from path, in INI style.")
      ("force-spike-list", po::value<std::string>(),
       "Read force spike list.")
	  ("synaptic-alpha-path", po::value<std::string>(),
	   "Set synaptic-alpha from path, only needed for DIF based model")
  ;
  // filter?

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    // Thanks stackoverflow
    auto basename = [](std::string filename) -> std::string {
      const size_t last_slash_idx = filename.find_last_of("\\/");
      if (std::string::npos != last_slash_idx) {
        filename.erase(0, last_slash_idx + 1);
      }
      return filename;
    };
    cout << "Usage: " << basename(argv[0]);
    cout << " --neuron-model arg";
    cout << " --net arg";
    cout << " [OPTION]..." << "\n";
    cout << "Neuron model simulator, with accurate spike timing and computation.\n\n";
    cout << desc << "\n";
    return 0;
  }

  if (vm.count("parameter-path")) {
    std::ifstream fin_opt(vm["parameter-path"].as<std::string>().c_str());
    po::store(po::parse_config_file(fin_opt, desc), vm);
  }
  po::notify(vm);

  // Set random seed for the global generator, mainly for Poisson event.
  if (vm.count("seed")) {
    const auto &v_str = vm["seed"].as< std::vector<std::string> >();
    std::vector<unsigned int> v;
    if (v_str[0] == "auto") {
      std::random_device rd;
      v.push_back(rd());
      v.push_back(rd());
    } else {
      v.resize(v_str.size());
      try {
        for (size_t k = 0; k < v_str.size(); k++) {
          v[k] = std::stoul(v_str[k]);
        }
      } catch (const std::invalid_argument e) {
        cerr << "Error: seed must be integer(s) or 'auto'.\n";
        return 1;
      }
    }
    std::seed_seq sseq(v.begin(), v.end());
    rand_eng.seed( sseq );
    if (vm.count("verbose-echo")) {
      cout << "Random seed:";
      for (auto val : v) cout << " " << val;
      cout << "\n";
    }
  }

  int rt = -1;
  try {
    rt = MainLoop(vm);
  } catch (const char *st) {
    cerr << "Error: " << st << "\n";
  }
//  fclose(fdbg);
  return rt;
}
