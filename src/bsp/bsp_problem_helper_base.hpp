#pragma once

#include "common.hpp"

namespace BSP {

  class BSPProblemHelperBase {
  public:
    virtual int get_state_dim() const = 0;
    virtual int get_control_dim() const = 0;
    virtual int get_observe_dim() const = 0;
    virtual int get_state_noise_dim() const = 0;
    virtual int get_observe_noise_dim() const = 0;
    virtual int get_belief_dim() const = 0;
    virtual int get_sigma_dof() const = 0;
    virtual int get_T() const = 0;

    BSPProblemHelperBase() {}
  };

  typedef boost::shared_ptr<BSPProblemHelperBase> BSPProblemHelperBasePtr;
}
