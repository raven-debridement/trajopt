#pragma once

#include "common.hpp"

namespace BSP {

  class BSPProblemHelper {
  public:
    int state_dim;
    int control_dim;
    int observe_dim;
    int state_noise_dim;
    int observe_noise_dim;
    int belief_dim;
    int sigma_dof;
  };

  typedef boost::shared_ptr<BSPProblemHelper> BSPProblemHelperPtr;
}
