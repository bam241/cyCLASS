#include "uncertainty.h"

namespace cyclass {

double get_corrected_param(double param, double param_uncertainty) {
  if (param_uncertainty == 0) {
    return param;
  } else {
    std::default_random_engine de(std::clock());
    std::normal_distribution<double> nd(param, param * param_uncertainty);
    double val_tmp = nd(de);
    double val = param;
    if ((double)val == (int)val_tmp) {
      val = std::round(val_tmp);
    } else {
      val = (double)val_tmp;
    }
    return val;
  }
}
double get_corrected_param(int param, double param_uncertainty) {
  if (param_uncertainty == 0) {
    return param;
  } else {
    std::default_random_engine de(std::clock());
    std::normal_distribution<double> nd(param, param * param_uncertainty);
    double val_tmp = nd(de);
    int val = param;
    if ((int)val == (int)val_tmp) {
      val = std::round(val_tmp);
    } else {
      val = (int)val_tmp;
    }
    return val;
  }
}
}
