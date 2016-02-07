#include "neuron_system_utils.h"

void FillNeuStateFromFile(TyNeuronalDymState &neu_dym_stat, const char *path)
{
  std::ifstream fin(path);
  double v;
  size_t j = 0;
  neu_dym_stat.dym_vals.setZero();
  while (true) {
    for (int k = 0; fin >> v, k < neu_dym_stat.dym_vals.cols(); k++) {
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
}

// Read network from text file. The number neurons should be known first.
void FillNetFromPath(TyNeuronalParams &pm, const std::string &name_net)
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
    // Read network from text file
    // Negative strength is possible, but that is unusual
    std::ifstream fin_net(name_net);
    double ev;
    for (int i = 0; i < n_neu; i++) {
      for (int j = 0; j < n_neu; j++) {
        fin_net >> ev;
        if (ev != 0 && std::isfinite(ev) && i != j) {
          net_coef.push_back(TyEdgeTriplet(i, j, ev));
        }
      }
    }
  }

  pm.net.setFromTriplets(net_coef.begin(), net_coef.end());
  for (int j = 0; j < pm.n_total(); j++) {
    if (pm.net.coeffRef(j,j)) {
      pm.net.coeffRef(j,j) = 0;
    }
  }
  pm.net.prune(std::numeric_limits<double>::min(), 1);  // remove zeros;
  pm.net.makeCompressed();
}

