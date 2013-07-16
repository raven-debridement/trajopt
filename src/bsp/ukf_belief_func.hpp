#pragma once

#include "common.hpp"
#include "traits.hpp"
#include "utils.hpp"
#include "state_func.hpp"
#include "observe_func.hpp"
#include "bsp_problem_helper_base.hpp"
#include "problem_state.hpp"
#include "types/belief_func.hpp"

namespace BSP {

  template< class _StateFuncT=StateFunc<>, class _ObserveFuncT=ObserveFunc<>, class _BeliefT=VectorXd >
  class UkfBeliefFunc : public BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> {
  public:
    define_belief_func_types();
    typedef boost::shared_ptr< UkfBeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT> > Ptr;

    UkfBeliefFunc() : BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>() {}

    UkfBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
      BeliefFunc<_StateFuncT, _ObserveFuncT, _BeliefT>(helper, f, h) {}

    virtual BeliefT operator()(const BeliefT& b, const ControlT& u, ObserveT* z_ptr=nullptr, ObserveT* observation_masks_ptr=nullptr, StateT* state_error_out=nullptr) const {
      // TODO
      
    }
  };
}
