#pragma once
#include "common.hpp"

namespace BSP {
  template< class _StateFuncT, class _ObserveFuncT, class _BeliefT>
  struct BeliefFuncTypes {
    typedef _StateFuncT StateFuncT;
    typedef _ObserveFuncT ObserveFuncT;
    typedef _BeliefT BeliefT;
    typedef typename StateFuncT::StateT StateT;
    typedef typename StateFuncT::ControlT ControlT;
    typedef typename StateFuncT::StateNoiseT StateNoiseT;
    typedef typename StateFuncT::StateGradT StateGradT;
    typedef typename StateFuncT::ControlGradT ControlGradT;
    typedef typename StateFuncT::StateNoiseGradT StateNoiseGradT;

    typedef typename ObserveFuncT::StateT ObserveStateT;
    typedef typename ObserveFuncT::ObserveT ObserveT;
    typedef typename ObserveFuncT::ObserveNoiseT ObserveNoiseT;
    typedef typename ObserveFuncT::ObserveStateGradT ObserveStateGradT;
    typedef typename ObserveFuncT::ObserveNoiseGradT ObserveNoiseGradT;

    typedef typename StateFuncT::Ptr StateFuncPtr;
    typedef typename ObserveFuncT::Ptr ObserveFuncPtr;

    typedef boost::shared_ptr< BeliefFunc<StateFuncT, ObserveFuncT, BeliefT> > Ptr;

    BOOST_STATIC_ASSERT( (boost::is_same< ObserveStateT, StateT >::value) );

    ENSURE_VECTOR(BeliefT);
    ENSURE_SAME_SCALAR_TYPE(StateT, ObserveT);
    ENSURE_SAME_SCALAR_TYPE(StateT, BeliefT);

    const static int _state_dim = StateFuncT::_state_dim;
    const static int _control_dim = StateFuncT::_control_dim;
    const static int _state_noise_dim = StateFuncT::_state_noise_dim;
    const static int _observe_dim = ObserveFuncT::_observe_dim;
    const static int _observe_noise_dim = ObserveFuncT::_observe_noise_dim;
    const static int _belief_dim = MatrixTraits<BeliefT>::rows;
    
    BOOST_STATIC_ASSERT( (_belief_dim >= _state_dim) );

    typedef typename StateFuncT::scalar_type scalar_type;

    typedef Matrix<scalar_type, _belief_dim, _belief_dim> BeliefGradT;
    typedef Matrix<scalar_type, _belief_dim, _control_dim> BeliefControlGradT;
    typedef Matrix<scalar_type, _state_dim, _state_dim> VarianceT;
    typedef Matrix<scalar_type, _observe_dim, _observe_dim> ObserveMatT;
    typedef Matrix<scalar_type, _state_dim, _observe_dim> KalmanT;
  };
}
