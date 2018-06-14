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
      //cerr << "Read " << j << " init data" << endl;
      break;
    }
  }
  if (j < neu_dym_stat.time_in_refractory.size()) {
    cerr << "No enough data read!" << endl;
    return -1;
  }
  return 0;
}

int FillVecFromFile(TyArrVals &v, const char *path)
{
  std::ifstream fin(path);
  if (fin.fail()) {
    cerr << "File not found: \"" << path << "\"" << endl;
    return -1;
  }
  v.clear();
  double d;
  while (fin >> d) {
    v.push_back(d);
  }
  return 0;
}

// Fill tau for conductance into dym_param from path.
// The order is tau_gE_s1, tau_gE, tau_gI_s1, tau_gI.
int FillTauG(TyNeuDymParam &dym_param, const char *path)
{
  TyArrVals v;
  int rt = FillVecFromFile(v, path);
  if (rt != 0) return rt;
  if (v.size() % 4 != 0) {
    cerr << "FillTauG: wrong content in file: \"" << path << "\"\n";
    return 1;
  }
  int n_neu = v.size() / 4;
  dym_param.resize(n_neu, 4);
  memcpy(dym_param.data(), v.data(), v.size()*sizeof(v[0]));
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
      char file_sig = fin_net.peek();
      if (file_sig == '#' || file_sig == '%') {
        std::string sig_line;
        std::getline(fin_net, sig_line);
        // remove sig line, assume correct sparse matrix.
      }
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


// YWS
// credit to http://www.cplusplus.com/forum/general/42594/
double string2double(const std::string& a)
{
	// Convert a string representation of a number into a floating point value.
	// Throws an int if the string contains anything but whitespace and a valid
	// numeric representation.
	//
	double result;
	std::string s(a);

	// Get rid of any trailing whitespace
	s.erase(s.find_last_not_of(" \f\n\r\t\v") + 1);

	// Read it into the target type
	std::istringstream ss(s);
	ss >> result;

	// Check to see that there is nothing left over
	if (!ss.eof())
		throw 1;

	return result;
}

bool isNumber(std::string str) {
	try
	{
		string2double(str);
	}
	catch (int)
	{
		return false;
	}
	return true;
}

void InitAlphaCoeffFromPath(TyNeuronalParams & pm, const std::string & name_coef, bool is_sparse)
{
	int n_neu = pm.n_total();

	if (isNumber(name_coef)) {
		double value = string2double(name_coef);
		TyArrVals onerow(n_neu, value);
		TyMatVals onemat(n_neu, onerow);

		pm.alpha.resize(n_neu, onemat);

	}
	else {
		std::ifstream fin_alpha(name_coef);
		if (!fin_alpha) {
			cerr << "Fail to open file! \"" << name_coef << "\"" << endl;
			throw "Fail to open file!\n";
		}
		if (!is_sparse) {
			// Read network from text file
			double a;
			for (int i = 0; i < n_neu; i++) {
				TyMatVals amat;

				for (int j = 0; j < n_neu; j++) {
					TyArrVals arow;
					for (int k = 0; k < n_neu; k++) {
						fin_alpha >> a;
						if (j == k && a != 0) {
							dbg_printf("Warning: in InitAlphaCoeffFromPath: it is unusual that alpha[i][i] to be nonzero. \n");
						}
						if (!std::isfinite(a)) {
							dbg_printf("Warning: in InitAlphaCoeffFromPath: it is unusual that alpha[i][i] to be infinite. \n");
						}
						arow.push_back(a);
					}
					amat.push_back(arow);
				}
				pm.alpha.push_back(amat);
			}

			if (!fin_alpha) {
				cerr << "Bad alpha coefficient file: \"" << name_coef << "\"\n"
					<< " Notice that alpha coefficent file shall contain n*n*n doubles (n is the number of neurons)\n";
				throw "Bad alpha coefficient file!\n";
			}
		}
		else {
			printf("DIF-GH model doesn't suppport sparse net");
			// TODO YWS

		}

	}
}


SparseMat ReadNetDelay(const std::string &dn_name, const SparseMat &net)
{
  std::ifstream fin_net(dn_name);
  int n_neu = net.outerSize();
  SparseMat dn(n_neu, n_neu);
  bool is_sparse = false;

  char file_sig = fin_net.peek();
  if (file_sig == '#' || file_sig == '%') {
    std::string sig_line;
    std::getline(fin_net, sig_line);
    std::transform(sig_line.begin(), sig_line.end(), sig_line.begin(), ::tolower);
    if (sig_line.find("sparse") != std::string::npos) {
      is_sparse = true;
    }
    size_t pos_size = sig_line.find("size");
    if (pos_size != std::string::npos) {
      auto fnum = [](char ch) -> bool { return '0'<=ch && ch<='9'; };
      auto frange = [&](std::string::const_iterator it0) {
        auto it1 = std::find_if(it0, sig_line.cend(), fnum);
        auto it2 = std::find_if_not(it1, sig_line.cend(), fnum);
        return std::make_pair(it1, it2);
      };
      auto rg1 = frange(sig_line.begin() + pos_size);
      if (rg1.first != sig_line.end()) {
        long n1 = std::stol(sig_line.substr(
              rg1.first-sig_line.begin(), rg1.second - rg1.first));
        auto rg2 = frange(rg1.second);
        long n2 = std::stol(sig_line.substr(
              rg2.first-sig_line.begin(), rg2.second - rg2.first));
        if (n1 != n_neu || n2 != n_neu) {
          cerr << "Number of neuron inconsistant:\n";
          cerr << "        net: "<< n_neu <<" x "<< n_neu <<"\n";
          cerr << "  delay-net: "<< n1 <<" x "<< n2 <<"\n";
        }
        cout << "  delay-net: "<< n1 <<" x "<< n2 <<"\n";
      }
    }
  }

  if (!is_sparse) {
    typedef Eigen::Triplet<double> TyEdgeTriplet;
    std::vector<TyEdgeTriplet> net_coef;

    // Read network delay from text file
    double ev;
    for (int i = 0; i < n_neu; i++) {
      for (int j = 0; j < n_neu; j++) {
        fin_net >> ev;
        // dn and net should have identical size
        if (net.coeff(i, j) != 0) {
          net_coef.push_back(TyEdgeTriplet(i, j, ev));
        }
      }
    }
    if (!fin_net) {
      cerr << "Bad network file: \"" << dn_name << "\"\n";
      throw "Bad network file!\n";
    }

    // Construct the sparse net delay
    dn.setFromTriplets(net_coef.begin(), net_coef.end());
    dn.makeCompressed();
  } else {
    dn = net;  // Make sure that the two networks have the same shape.
    for (int k=0; k<dn.outerSize(); k++)
      for (SparseMat::InnerIterator it(dn, k); it; ++it) {
        it.valueRef() = 0;
      }

    int i, j;
    double ev;
    while (fin_net >> i >> j >> ev) {
      if (ev != 0 && std::isfinite(ev) && i != j) {
        if (i > n_neu || j > n_neu) {
          cout << "i j ev = " << i << " " << j << " " << ev << "\n";
          cout << "n_neuron = " << n_neu << "\n";
          throw "Number of neurons inconsistant.";
        }
        if (i <= 0 || j <= 0) {
          cout << "i j ev = " << i << " " << j << " " << ev << "\n";
          throw "Neuron index start from 1.";
        }
        if (net.coeff(i-1, j-1) != 0) {
          dn.coeffRef(i-1, j-1) = ev;
        }
      }
    }
    dn.makeCompressed();
    if (! fin_net.eof() && (fin_net.rdstate() & std::ios_base::failbit)) {
      cerr << "Bad network file: \"" << dn_name << "\"\n";
      throw "Bad network file!\n";
    }
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
