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

#define NDEBUG  // disable assert() and disable checks in Eigen

#include <cassert>
#include "common_header.h"
#include "legacy_lib.h"

#include "single_neuron_dynamics.h"
#include "neuron_population.h"
#include "simulator_exact_order.h"

#include "neuron_population_cont_synaptic.h"
#include "simulator_cont_synaptic.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef std::vector<unsigned int> VecUInt;

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
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
    cerr << "Error: Neuron model not specified. See --help.\n";
    return -1;
  }
  const std::string &str_nm = vm["neuron-model"].as<std::string>();

  bool b_verbose = false;
  if (vm.count("verbose")) {
    b_verbose = true;
  }

  // Neuron models.
  enum EnumNeuronModel {
    LIF_G, LIF_GH, HH_G, HH_GH, HH_PT_GH, HH_FT_GH,
    HH_G_sine, HH_GH_sine, HH_PT_GH_sine, HH_FT_GH_sine,
    HH_G_extI, HH_GH_extI,
    HH_GH_cont_syn
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
  }  else if (str_nm ==    "HH-GH-cont-syn") {
    enum_neuron_model =     HH_GH_cont_syn;
  } else {
    cerr << "Unrecognized neuron model. See --help.\n";
    return -1;
  }

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

  if (vm.count("pr-mul")) {
    const auto &v = vm["pr-mul"].as< std::vector<double> >();
    if (v.size() != (size_t)n_neu) {
      cerr << "parameter --pr-mul should have the same number of terms as neuron number.\n";
      return -1;
    }
    for (int i = 0; i < n_neu; i++) {
      pm.arr_pr[i] *= v[i];
    }
  }

  if (vm.count("ps-mul")) {
    const auto &v = vm["ps-mul"].as< std::vector<double> >();
    if (v.size() != (size_t)n_neu) {
      cerr << "parameter --ps-mul should have the same number of terms as neuron number.\n";
      return -1;
    }
    for (int i = 0; i < n_neu; i++) {
      pm.arr_ps[i] *= v[i];
    }
  }

  if (vm.count("psi") || vm.count("pri")) {
    cerr << "option --psi and --pri are not supported yet!" << endl;
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
    return -1;
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
      default:
        cerr << "This delay net type is not supported yet.\n";
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
    }
    p_neu_pop->SetSynapticDelay(vm["synaptic-delay"].as<double>());
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
  if (str_simu_method.size() > 0 && str_simu_method != "auto") {
    str_simu_method = vm["simulation-method"].as<std::string>();
    if (str_simu_method == "simple") {
      p_neu_simu = new NeuronSimulatorSimple(pm, e_dt);
    } else if (str_simu_method == "SSC") {  // Spike-Spike-Correction
      p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt);
    } else if (str_simu_method == "SSC-Sparse") {
      p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse(pm, e_dt);
    } else if (str_simu_method == "SSC-Sparse2") {
      p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse2(pm, e_dt);
    } else if (str_simu_method == "big-net-delay") {
      p_neu_simu = new NeuronSimulatorBigNetDelay(pm, e_dt);
    } else if (str_simu_method == "big-delay") {
      p_neu_simu = new NeuronSimulatorBigDelay(pm, e_dt);
    } else if (str_simu_method == "cont-syn") {
      p_neu_simu = new NeuronSimulatorCont(pm, e_dt);
    } else {
      cerr << "No this simulation method:\"" << str_simu_method << "\"\n";
      return -2;
    }
  } else {
    // Default simulator
    if (HH_GH_cont_syn == enum_neuron_model) {
      str_simu_method = "cont-syn";
      p_neu_simu = new NeuronSimulatorCont(pm, e_dt);
    } else if (vm.count("synaptic-net-delay")) {
      str_simu_method = "big-net-delay";
      p_neu_simu = new NeuronSimulatorBigNetDelay(pm, e_dt);
    } else if (vm.count("synaptic-delay")) {
      str_simu_method = "big-delay";
      p_neu_simu = new NeuronSimulatorBigDelay(pm, e_dt);
    } else {
      str_simu_method = "SSC";
      p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt);
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
  
  if (vm.count("no-threshold")) {
    p_neu_pop->DisableThreshold();
  }
  
  bool b_init_loaded = false;
  if (vm.count("initial-state-path")) {
    int rt = FillNeuStateFromFile(p_neu_pop->GetDymState(),
                 vm["initial-state-path"].as<std::string>().c_str());
    if (rt == 0) {
      b_init_loaded = true;
      cout << "Initial state loaded." << endl;
    }
    //p_neu_simu->SaneTestVolt();
  }
  if (!b_init_loaded) {
    // Fill with default values
    const double *v = p_neuron_model->Get_dym_default_val();
    auto &ds = p_neu_pop->GetDymState().dym_vals;
    int n_var = p_neuron_model->Get_n_dym_vars();
    for (int i = 0; i < p_neu_pop->n_neurons(); i++) {
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
                              vm["input-event-path"].as<std::string>().c_str());
    //cout << "input event loaded!" << endl;  // Debug info
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
  
  // 
  TySpikeEventVec force_spike_list;
  if (vm.count("force-spike-list")) {
    ReadSpikeList(force_spike_list, 
      vm["force-spike-list"].as<std::string>().c_str());
    printf("ReadSpikeList l=%lu\n", force_spike_list.size());
  }

  // Function for save data to file
  auto func_save_dym_values = [
    b_output_volt, &fout_volt,
    b_output_G, &fout_G,
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
    if (b_output_G) {
      int id_gE = p_neuron_model->Get_id_gE();
      int id_gI = p_neuron_model->Get_id_gI();
      std::vector<double> v(2*neu_pop.n_neurons());
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        v[2*j  ] = neu_pop.GetDymState().dym_vals(j, id_gE);
        v[2*j+1] = neu_pop.GetDymState().dym_vals(j, id_gI);
      }
      fout_G.write((char*)v.data(), 2*neu_pop.n_neurons()*sizeof(double));
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

  std::vector<size_t> vec_n_spike(p_neu_pop->n_neurons());  // count the number of spikes
  TySpikeEventVec ras;                            // record spike raster
  int n_dt_in_stv = int(e_stv / e_dt + 0.1);
  int count_n_dt_in_stv = n_dt_in_stv;
  size_t n_step = (size_t)(e_t / e_dt);

  // Show progress
  if (b_verbose) {
    printf("  Initialization     : %3.3f s\n", toc(t_begin));
    fflush(stdout);
    printf("  Simulation         : ");
  }

  // Main loop
  t_begin = tic();
  int progress_percent = 0;
  int str_len = 0;
  for (size_t i = 0; i < n_step; i++) {
    ForceSpikeOnList(p_neu_pop, force_spike_list, p_neu_simu->GetT(), e_dt);
    p_neu_simu->NextDt(p_neu_pop, ras, vec_n_spike);
    
    //p_neu_simu->SaneTestState();
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
    func_save_dym_values(*p_neu_pop);
    
    // Show progress
    if (b_verbose && (i+1.0)/n_step >= (progress_percent+1.0)/100) {
      progress_percent = (int)((i+1.0)/n_step * 100);
      // For matlab compatible (otherwise just printf('\r') will be fine)
      for (int j=0; j<str_len; j++) printf("\b");
      str_len = printf("%3.3f s (%d%%)",
                       toc(t_begin), progress_percent);
      fflush(stdout);
    }
  }
  if (b_verbose) {
    printf("\n");
    fflush(stdout);
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
  // Declare the supported options.
  po::options_description desc("Options");
  // http://stackoverflow.com/questions/3621181/short-options-only-in-boostprogram-options
  desc.add_options()
      ("neuron-model",  po::value<std::string>(),
       "One of LIF-G, LIF-GH, HH-G, HH-GH, HH-PT-GH, HH-FT-GH, HH-G-sine, HH-GH-sine, HH-PT-GH-sine, HH-FT-GH-sine, HH-G-extI, HH-GH-extI, HH-GH-cont-syn.")
      ("simulation-method",  po::value<std::string>(),
       "One of simple, SSC, SSC-Sparse, SSC-Sparse2, big-delay, big-net-delay, cont-syn, auto. Some combinations of neuron model and simulator are mutually exclusive, hence not allowed. If not specify, a suitable simulator will be choosen automatically.")
      ("help,h",
       "Produce help message.")
      ("verbose,v",
       "Show progress.")
      ("t",    po::value<double>()->default_value(1e3),
       "Simulation time, in ms.")
      ("dt",   po::value<double>()->default_value(1.0/32),
       "Simulation delta t (dt, time step), in ms.")
      ("stv",  po::value<double>(),
       "Output sampling interval. Must be multiples of dt. Default set to dt.")
      ("nE",   po::value<unsigned int>()->default_value(1),
       "Number of excitatory neurons.")
      ("nI",   po::value<unsigned int>()->default_value(0),
       "Number of inhibitory neurons.")
      ("net",  po::value<std::string>(),
       "Network file path.")
      ("scee", po::value<double>()->default_value(0.0),
       "Cortical strength E->E.")
      ("scie", po::value<double>()->default_value(0.0),
       "Cortical strength E->I.")
      ("scei", po::value<double>()->default_value(0.0),
       "Cortical strength I->E.")
      ("scii", po::value<double>()->default_value(0.0),
       "Cortical strength I->I.")
      ("ps",   po::value<double>(),
       "Poisson input strength.")
      ("pr",   po::value<double>()->default_value(1.0),
       "Poisson input rate, in 1/ms.")
      ("psi",  po::value<double>(),
       "Poisson input strength, inhibitory type.")
      ("pri",  po::value<double>(),
       "Poisson input rate, inhibitory type.")
      ("pr-mul", po::value< std::vector<double> >()->multitoken(),
       "Poisson input rate multiper.")
      ("ps-mul", po::value< std::vector<double> >()->multitoken(),
       "Poisson input strength multiper.")
      ("seed",   po::value< VecUInt >()->multitoken()->default_value(VecUInt{1}, "1"),
       "Random seed for Poisson events. One or several unsigned integers (0 ~ 2^32-1).")
      ("seed-auto",
       "Auto set random seed. This option overrides --seed.")
      ("no-threshold",
       "Disable the threshold detection and voltage reset. Note that HH based neurons will still fire.")
      ("current-sine-amp",         po::value<double>(),
       "Set the current input sine amplitude.")
      ("current-sine-freq", po::value<double>(),
       "Set the current input sine frequency, in kHz.")
      ("extI", po::value<double>(),
       "Set the custom current input extI.")
      ("synaptic-delay", po::value<double>(),
       "Set a synaptic delay for the network.")
      ("synaptic-net-delay", po::value<std::string>(),
       "Set the synaptic delays for the network. Path to the matrix text file.")
      ("volt-path,o",      po::value<std::string>(),
       "Output volt to path. In raw binary format.")
      ("ras-path",         po::value<std::string>(),
       "Output spike events to path.")
      ("isi-path",         po::value<std::string>(),
       "Output mean InterSpike Interval to path.")
      ("conductance-path", po::value<std::string>(),
       "Output conductance to path.")
      ("ion-gate-path", po::value<std::string>(),
       "Output gating variables to path. In raw binary format.")
      ("output-first-data-point",
       "Also output initial condition in volt-path and conductance-path.")
      ("initial-state-path", po::value<std::string>(),
       "Read initial state from path.")
      ("input-event-path", po::value<std::string>(),
       "Read input event from path.")
      ("parameter-path", po::value<std::string>(),
       "Read parameters from path, in INI style.")
      ("force-spike-list", po::value<std::string>(),
       "Read force spike list.")
  ;
  // ps-mul
  // verbose : for progress percentage
  // filter
  // vector input

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
    cout << " --neuron-model [arg]";
    cout << " --net [arg]";
    cout << " --ps [arg]";
    cout << " [OPTION]..." << "\n";
    cout << "Neuron model simulator, with accurate firing timing and computation.\n";
    cout << desc << "\n";
    return 1;
  }

  if (vm.count("parameter-path")) {
    std::ifstream fin_opt(vm["parameter-path"].as<std::string>().c_str());
    po::store(po::parse_config_file(fin_opt, desc), vm);
  }
  po::notify(vm);

  // Set random seed for the global generator, mainly for Poisson event.
  if (vm.count("seed")) {
    const auto &v = vm["seed"].as< VecUInt >();
//    for (size_t k = 0; k < v.size(); k++)
//      cout << "v[" << k <<"] = " << v[k] << endl;
    std::seed_seq sseq(v.begin(), v.end());
    rand_eng.seed( sseq );
  }
  if (vm.count("seed-auto")) {
    std::random_device rd;
    std::seed_seq sseq{rd(), rd()};
    rand_eng.seed( sseq );
  }

  int rt = -1;
  try {
    rt = MainLoop(vm);
  } catch (const char *st) {
    cerr << "Error: " << st << "\n";
  }

  return rt;
}
