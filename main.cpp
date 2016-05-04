// g++ -g -O0 -std=c++11 -Wall main.cpp math_helper.cpp -lboost_program_options -o bin/vec_IFsimu
//Could faster:
// g++ -g -O2 -falign-functions=16 -falign-loops=16

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
   * possibility to add synaptic delay? add to Poisson event queue?
   * isomerisom of neurons in a network?
*/

#define NDEBUG  // disable assert() and disable checks in Eigen

#include <libgen.h>  // for basename
#include <cassert>
#include "common_header.h"
#include "legacy_lib.h"

#include "single_neuron_dynamics.h"
//#include "neuron_system_utils.h"
#include "neuron_population.h"
#include "simulator_exact_order.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

int MainLoop(const po::variables_map &vm)
{
  if (!vm.count("neuron-model")) {
    cerr << "Error: Neuron model not specified. See --help.\n";
    return -1;
  }
  Ty_Neuron_Dym_Base *p_neuron_model;
  const std::string &str_nm = vm["neuron-model"].as<std::string>();
  enum EnumNeuronModel {LIF_G, LIF_GH, HH_GH, HH_GH_sine };
  EnumNeuronModel enum_neuron_model;
  if (str_nm == "LIF-G" || str_nm == "LIF-G-Sparse") {
    p_neuron_model = new Ty_LIF_G();
    enum_neuron_model = LIF_G;
  } else if (str_nm == "LIF-GH" || str_nm == "LIF-GH-Sparse") {
    p_neuron_model = new Ty_LIF_GH();
    enum_neuron_model = LIF_GH;
  } else if (str_nm == "HH-GH" || str_nm == "HH-GH-Sparse") {
    p_neuron_model = new Ty_HH_GH();
    enum_neuron_model = HH_GH;
  }  else if (str_nm == "HH-GH-sine" || str_nm == "HH-GH-sine-Sparse") {
    p_neuron_model = new Ty_HH_GH_sine();
    enum_neuron_model = HH_GH_sine;
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

  NeuronPopulationBase * p_neu_pop = nullptr;
  switch (enum_neuron_model) {
    case LIF_G:
      p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_LIF_G>(pm);
      break;
    case LIF_GH:
      p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_LIF_GH>(pm);
      break;
    case HH_GH:
      p_neu_pop = new NeuronPopulationDeltaInteractTemplate<Ty_HH_GH>(pm);
      break;
    case HH_GH_sine:
      p_neu_pop = new NeuronPopulationDeltaInteractSine<Ty_HH_GH_sine>(pm);
      break;
  }
  //NeuronPopulationDeltaInteract neu_pop(p_neuron_model, pm);

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
  
  // Parameters about current input
  if (enum_neuron_model == HH_GH_sine) {
    auto p_neu_pop_sine = static_cast<
      NeuronPopulationDeltaInteractSine<Ty_HH_GH_sine> *>(p_neu_pop);
    if (vm.count("sine-current-amplitude")) {
      p_neu_pop_sine->SetSineAmplitude(
          vm["sine-current-amplitude"].as<double>());
    }
    if (vm.count("sine-current-angular-frequency")) {
      p_neu_pop_sine->SetSineAngularFrequency(
          vm["sine-current-angular-frequency"].as<double>());
    }
  }

  auto fout_try_open = [&vm](const char * const st_id, std::ofstream &fs)
    -> bool {
      if (vm.count(st_id)) {
        std::string st_path = vm[st_id].as<std::string>();
        CheckDirAndCreate( st_path );
        fs.open( st_path );
        if (!fs) {
          cerr << "Error: Failed to open file \"" << st_path << "\" for output." << "\n";
          throw "Failed to open file";
        }
        return true;
      } else {
        return false;
      }
    };

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

  // simulator for the neural network
  NeuronSimulatorBase *p_neu_simu = nullptr;
  p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt);
  /*
  if (str_nm == "LIF-G" || str_nm == "LIF-GH" || str_nm == "HH-GH") {
    p_neu_simu = new NeuronSimulatorExactSpikeOrder(pm, e_dt);
  } else if (str_nm == "LIF-G-Sparse" || str_nm == "LIF-GH-Sparse"
             || str_nm == "HH-GH-Sparse") {
    p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse(pm, e_dt);
  } else if (str_nm == "LIF-G-Sparse2" || str_nm == "LIF-GH-Sparse2"
             || str_nm == "HH-GH-Sparse2") {
    p_neu_simu = new NeuronSimulatorExactSpikeOrderSparse2(pm, e_dt);
  }*/

  if (vm.count("initial-state-path")) {
    FillNeuStateFromFile(p_neu_pop->GetDymState(),
                         vm["initial-state-path"].as<std::string>().c_str());
    cout << "initial state loaded!" << endl;
    p_neu_simu->SaneTestVolt();
  }

  if (vm.count("input-event-path")) {
    FillPoissonEventsFromFile(p_neu_simu->Get_poisson_time_vec(),
                              vm["input-event-path"].as<std::string>().c_str());
    cout << "input event loaded!" << endl;
  }

  auto func_save_dym_values = [
    b_output_volt, &fout_volt,
    b_output_G, &fout_G,
    b_output_HH_gate, &fout_HH_gate,
    p_neuron_model]
    (const NeuronPopulationBase &neu_pop) {
    if (b_output_volt) {
      int id_V = p_neuron_model->Get_id_V();
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        double d = neu_pop.GetDymState().dym_vals(j, id_V);
        fout_volt.write((char*)&d, sizeof(double));
      }
    }
    if (b_output_G) {
      int id_gE = p_neuron_model->Get_id_gE();
      int id_gI = p_neuron_model->Get_id_gI();
      for (int j = 0; j < neu_pop.n_neurons(); j++) {
        double d1 = neu_pop.GetDymState().dym_vals(j, id_gE);
        double d2 = neu_pop.GetDymState().dym_vals(j, id_gI);
        fout_G.write((char*)&d1, sizeof(double));
        fout_G.write((char*)&d2, sizeof(double));
      }
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

  // Main loop
  for (size_t i = 0; i < n_step; i++) {
    p_neu_simu->NextDt(p_neu_pop, ras, vec_n_spike);
    p_neu_simu->SaneTestState();
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
  }

  if (vm.count("isi-path")) {
    std::ofstream fout_isi( vm["isi-path"].as<std::string>() );
    fout_isi.precision(17);
    for (int j = 0; j < p_neu_pop->n_neurons(); j++) {
      fout_isi << e_t / vec_n_spike[j] << '\t';
    }
  }

  return 0;
}

int main(int argc, char *argv[])
{
  // Declare the supported options.
  po::options_description desc("Options");
  // http://stackoverflow.com/questions/3621181/short-options-only-in-boostprogram-options
  desc.add_options()
      ("neuron-model",  po::value<std::string>(),
       "One of LIF-G, LIF-GH, HH-GH and HH-GH-sine.")
      ("help,h",
       "Produce help message.")
      ("verbose,v",
       "Show progress.")
      ("t",    po::value<double>()->default_value(1e3),
       "Simulation time, in ms.")
      ("dt",   po::value<double>()->default_value(1.0/32),
       "Simulation delta t (dt, time step), in ms.")
      ("stv",  po::value<double>(),
       "Delta t for output (sampling). Must be multiples of dt.")
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
       "Poisson input rate, unit: 1/ms.")
      ("psi",  po::value<double>(),
       "Poisson input strength, inhibitory type.")
      ("pri",  po::value<double>(),
       "Poisson input rate, inhibitory type.")
      ("pr-mul", po::value<double>()->multitoken(),
       "Poisson input rate multiper.")
      ("seed",   po::value<unsigned int>()->default_value(1),
       "Random seed for Poisson events. An unsigned integer (0~2^32-1).")
      ("seed-auto",
       "Auto set random seed. This option overrides --seed.")
      ("sine-current-amplitude",         po::value<double>(),
       "Set the sine amplitude for current input.")
      ("sine-current-angular-frequency", po::value<double>(),
       "Set the sine angular frequency for current input.")
      ("volt-path,o",      po::value<std::string>(),
       "Volt output file path.")
      ("ras-path",         po::value<std::string>(),
       "RAS output file path.")
      ("isi-path",         po::value<std::string>(),
       "ISI output file path.")
      ("conductance-path", po::value<std::string>(),
       "Conductance output file path.")
      ("ion-gate-path", po::value<std::string>(),
       "Gating variables output file path.")
      ("initial-state-path", po::value<std::string>(),
       "Initial state file path.")
      ("input-event-path", po::value<std::string>(),
       "Input event file path.")
      ("output-first-data-point",
       "Also output initial condition in volt-path and conductance-path.")
      ("parameter-path", po::value<std::string>(),
       "Path for parameters, in INI style.")
  ;
  // ps-mul
  // verbose : for progress percentage
  // filter
  // vector input

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    cout << "Usage: " << basename(argv[0]) << " [OPTION]..." << "\n";
    cout << "Neuron model simulator, with accurate firing timing and computation.\n";
    cout << desc << "\n";
    return 1;
  }

  if (vm.count("parameter-path")) {
    std::ifstream fin_opt(vm["parameter-path"].as<std::string>().c_str());
    po::store(po::parse_config_file(fin_opt, desc), vm);
  }
  po::notify(vm);

  // Set random seed for Poisson event generator.
  if (vm.count("seed")) {
    rand_eng.seed( vm["seed"].as<unsigned int>() );
  }
  if (vm.count("seed-auto")) {
    std::random_device rd;
    std::seed_seq sseq{rd(), rd()};
    rand_eng.seed( sseq );
  }

  // Select neuron model.
  int rt = MainLoop(vm);

  return rt;
}
