#include "common.hpp"
#include "utils.hpp"
#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
  #include <direct.h>
  #define GetCurrentDir _getcwd
#else
  #include <unistd.h>
  #define GetCurrentDir getcwd
#endif

namespace BSP {
  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    int n_arr = name_prefix.size();
    assert(n_arr == newvars.size());

    vector<MatrixXi> index(n_arr);
    for (int i=0; i < n_arr; ++i) {
      newvars[i]->resize(rows, cols[i]);
      index[i].resize(rows, cols[i]);
    }

    vector<string> names;
    vector<double> all_lbs;
    vector<double> all_ubs;
    int var_idx = prob.getNumVars();
    for (int i=0; i < rows; ++i) {
      for (int k=0; k < n_arr; ++k) {
        for (int j=0; j < cols[k]; ++j) {
          index[k](i,j) = var_idx;
          names.push_back( (boost::format("%s_%i_%i")%name_prefix[k]%i%j).str() );
          all_lbs.push_back(lbs[k]);
          all_ubs.push_back(ubs[k]);
          ++var_idx;
        }
      }
    }
    prob.createVariables(names, all_lbs, all_ubs); // note that w,r, are both unbounded

    const vector<Var>& vars = prob.getVars();
    for (int k=0; k < n_arr; ++k) {
      for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols[k]; ++j) {
          (*newvars[k])(i,j) = vars[index[k](i,j)];
        }
      }
    }
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    vector<double> lbs(newvars.size(), -INFINITY);
    vector<double> ubs(newvars.size(), INFINITY);
    AddVarArrays(prob, rows, cols, lbs, ubs, name_prefix, newvars);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars) {
    vector<VarArray*> arrs(1, &newvars);
    vector<string> prefixes(1, name_prefix);
    vector<int> colss(1, cols);
    vector<double> lbs(1, lb);
    vector<double> ubs(1, ub);
    AddVarArrays(prob, rows, colss, lbs, ubs, prefixes, arrs);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars) {
    AddVarArray(prob, rows, cols, -INFINITY, INFINITY, name_prefix, newvars);
  }

  void seed_random() {
    unsigned int randval;
    FILE *f;
    f = fopen("/dev/random", "r");
    fread(&randval, sizeof(randval), 1, f);
    fclose(f);
    srand(randval);
  }

  double sigmoid(double x) {
    return 1. / (1. + exp(-x));
  }

  string get_current_directory() {
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
      throw std::runtime_error("cannot get current path");
    }
    return string(cCurrentPath);
  }

}
