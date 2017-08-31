#include "neuron_system_utils.h"

int FillNeuStateFromFile(TyNeuronalDymState &neu_dym_stat, const char *path)
{
  std::ifstream fin(path);
  if (fin.fail()) {
    cerr << "Initial value file not found! \"" << path << "\"" << endl;
    return -1;
  }

  double v;
  size_t j = 0;
  neu_dym_stat.dym_vals.setZero();
  while (true) {
    for (int k = 0; k < neu_dym_stat.dym_vals.cols(); k++) {
      fin >> v;
      neu_dym_stat.dym_vals(j, k) = v;
    }
    if (fin.fail())
      break;
    neu_dym_stat.time_in_refractory[j] = 0;
    j++;
    if (j >= neu_dym_stat.time_in_refractory.size()) {
      cerr << "read " << j << " init data" << endl;
      break;
    }
  }
  if (j < neu_dym_stat.time_in_refractory.size()) {
    cerr << "No enough data read!" << endl;
    return -1;
  }
  return 0;
}

// Read network from text file. The number neurons should be known first.
void FillNetFromPath(TyNeuronalParams &pm, const std::string &name_net,
                     bool is_sparse)
{
  typedef Eigen::Triplet<double> TyEdgeTriplet;
  std::vector<TyEdgeTriplet> net_coef;
  int n_neu = pm.n_total();

  if (name_net == "-") {
    // Fully connected network
    for (int i = 0; i < n_neu; i++) {
      for (int j = 0; j < n_neu; j++) {
        if (i==j) {
          continue;
        }
        net_coef.push_back(TyEdgeTriplet(i, j, 1.0));
      }
    }
  } else if (name_net == "--") {
    // Simple ring structure
    net_coef.push_back(TyEdgeTriplet(0, n_neu-1, 1.0));
    for (int i = 1; i < n_neu; i++) {
      net_coef.push_back(TyEdgeTriplet(i, i-1, 1.0));
    }
  } else {
    std::ifstream fin_net(name_net);
    if (!fin_net) {
      cerr << "Fail to open file! \"" << name_net << "\"" << endl;
      throw "Fail to open file!\n";
    }
    if (!is_sparse) {
      // Read network from text file
      // Negative strength is possible, but that is unusual
      double ev;
      for (int i = 0; i < n_neu; i++) {
        for (int j = 0; j < n_neu; j++) {
          fin_net >> ev;
          if (ev != 0 && std::isfinite(ev) && i != j) {
            net_coef.push_back(TyEdgeTriplet(i, j, ev));
          }
        }
      }
      if (!fin_net) {
        cerr << "Bad network file: \"" << name_net << "\"\n";
        throw "Bad network file!\n";
      }
    } else {
      int i, j;
      double ev;
      while (fin_net >> i >> j >> ev) {
        if (ev != 0 && std::isfinite(ev) && i != j) {
          if (i > n_neu || j > n_neu) {
            cout << "i j ev = " << i << " " << j << " " << ev << "\n";
            cout << "n_neuron = " << n_neu << "\n";
            throw "Number of neurons inconsistant.";
          }
          net_coef.push_back(TyEdgeTriplet(i-1, j-1, ev));
        }
      }
      if (! fin_net.eof() && (fin_net.rdstate() & std::ios_base::failbit)) {
        cerr << "Bad network file: \"" << name_net << "\"\n";
        throw "Bad network file!\n";
      }
    }
  }
  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());
  for (int j = 0; j < pm.n_total(); j++) {
    if (pm.net.coeff(j,j)) {
      pm.net.coeffRef(j,j) = 0;
    }
  }
  pm.net.prune(std::numeric_limits<double>::min(), 1);  // remove zeros;
  pm.net.makeCompressed();

  if (pm.net.outerSize() != n_neu) {
    cerr << "Error: number of neurons inconsistant (network file and command line)!\n";
    throw "Number of neurons inconsistant";
  }
}

SparseMat ReadNetDelay(const std::string &dn_name, const SparseMat &net)
{
  typedef Eigen::Triplet<double> TyEdgeTriplet;
  std::vector<TyEdgeTriplet> net_coef;

  // Read network delay from text file
  std::ifstream fin_net(dn_name);
  double ev;
  int n_neu = net.outerSize();
  for (int i = 0; i < n_neu; i++) {
    for (int j = 0; j < n_neu; j++) {
      fin_net >> ev;
      // dn and net should have identical size
      if (net.coeff(i, j) != 0) {
        net_coef.push_back(TyEdgeTriplet(i, j, ev));
      }
    }
  }

  // Construct the sparse net delay
  SparseMat dn(n_neu, n_neu);
  dn.setFromTriplets(net_coef.begin(), net_coef.end());
  dn.makeCompressed();

  if (!fin_net) {
    cerr << "Bad network delay file was read!\n";
    throw "Bad network file was read!\n";
  }

  return dn;
}

int ReadSpikeList(TySpikeEventVec &spike_list, const char *path)
{
  std::ifstream fin(path);
  if (fin.fail()) {
    cerr << "Fail to open file! \"" << path << "\"" << endl;
    return -1;
  }
  spike_list.resize(0);
  size_t id;
  double time;
  while (fin >> id >> time) {
    if (id == 0) {
      cerr << "ReadSpikeList(): Read id=0, neuron id start from 1.\n";
      return -1;
    }
    spike_list.emplace_back(time, id - 1);
  }
  return 0;
}
