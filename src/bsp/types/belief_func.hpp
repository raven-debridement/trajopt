#pragma once

#define define_belief_func_types()                                                                                                    \
  typedef _StateFuncT StateFuncT;                                                                                                     \
  typedef _ObserveFuncT ObserveFuncT;                                                                                                 \
  typedef _BeliefT BeliefT;                                                                                                           \
  typedef typename StateFuncT::StateT StateT;                                                                                         \
  typedef typename StateFuncT::ControlT ControlT;                                                                                     \
  typedef typename StateFuncT::StateNoiseT StateNoiseT;                                                                               \
  typedef typename StateFuncT::StateGradT StateGradT;                                                                                 \
  typedef typename StateFuncT::ControlGradT ControlGradT;                                                                             \
  typedef typename StateFuncT::StateNoiseGradT StateNoiseGradT;                                                                       \
                                                                                                                                      \
  typedef typename ObserveFuncT::StateT ObserveStateT;                                                                                \
  typedef typename ObserveFuncT::ObserveT ObserveT;                                                                                   \
  typedef typename ObserveFuncT::ObserveNoiseT ObserveNoiseT;                                                                         \
  typedef typename ObserveFuncT::ObserveStateGradT ObserveStateGradT;                                                                 \
  typedef typename ObserveFuncT::ObserveNoiseGradT ObserveNoiseGradT;                                                                 \
                                                                                                                                      \
  typedef typename StateFuncT::Ptr StateFuncPtr;                                                                                      \
  typedef typename ObserveFuncT::Ptr ObserveFuncPtr;                                                                                  \
                                                                                                                                      \
  static_assert(std::is_same<StateT, ObserveStateT>::value, "observe function and state function must have the same state type");     \
  static_assert(is_vector<BeliefT>::value, "belief must be an Eigen3 vector");                                                        \
                                                                                                                                      \
  static_assert(is_same_scalar_type<StateT, ObserveT>::value, "state and observe vectors must have the same scalar type");            \
  static_assert(is_same_scalar_type<StateT, BeliefT>::value, "state and belief vectors must have the same scalar type");              \
                                                                                                                                      \
  const static int _state_dim = StateFuncT::_state_dim;                                                                               \
  const static int _control_dim = StateFuncT::_control_dim;                                                                           \
  const static int _state_noise_dim = StateFuncT::_state_noise_dim;                                                                   \
  const static int _observe_dim = ObserveFuncT::_observe_dim;                                                                         \
  const static int _observe_noise_dim = ObserveFuncT::_observe_noise_dim;                                                             \
  const static int _belief_dim = MatrixTraits<BeliefT>::rows;                                                                         \
                                                                                                                                      \
  static_assert( (_belief_dim == -1 || _belief_dim >= _state_dim), "belief vector should be at least as large as the state vector" ); \
                                                                                                                                      \
  typedef typename StateFuncT::scalar_type scalar_type;                                                                               \
                                                                                                                                      \
  typedef Matrix<scalar_type, _belief_dim, _belief_dim> BeliefGradT;                                                                  \
  typedef Matrix<scalar_type, _state_dim, _belief_dim> SigmaPointsGradT;                                                                  \
  typedef Matrix<scalar_type, _belief_dim, _control_dim> BeliefControlGradT;                                                          \
  typedef Matrix<scalar_type, _state_dim, _state_dim> VarianceT;                                                                      \
  typedef Matrix<scalar_type, _observe_dim, _observe_dim> ObserveMatT;                                                                \
  typedef Matrix<scalar_type, _state_dim, _observe_dim> KalmanT;
