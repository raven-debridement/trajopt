#pragma once

#include "common.hpp"
#include "utils/logging.hpp"
#include "sco/optimizers.hpp"

namespace BSP {

  class BSPTrustRegionSQP : public BasicTrustRegionSQP {
  public:

    BSPTrustRegionSQP();
    BSPTrustRegionSQP(OptProbPtr prob);
    void addMeritDoneCallback(const Callback& cb);
    OptStatus optimize();
  protected:
    vector<Callback> merit_done_callbacks_;
    void callMeritDoneCallbacks(DblVec& x);
  };

}
