#pragma once

namespace BSP {
  class ProblemState {
  public:
    int state_dim;
    int control_dim;
    int observe_dim;
    int state_noise_dim;
    int observe_noise_dim;
    int belief_dim;
    int sigma_dof;
    ProblemState() {}
    ProblemState(BSPProblemHelperBasePtr helper) {
      state_dim = helper->get_state_dim();
      control_dim = helper->get_control_dim();
      observe_dim = helper->get_observe_dim();
      state_noise_dim = helper->get_state_noise_dim();
      observe_noise_dim = helper->get_observe_noise_dim();
      belief_dim = helper->get_belief_dim(); 
      sigma_dof = helper->get_sigma_dof(); 
    }
  };
}
