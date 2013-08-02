#pragma once

#define define_observe_func_types()                                                                                            \
  typedef _StateT StateT;                                                                                                      \
  typedef _ObserveT ObserveT;                                                                                                  \
  typedef _ObserveNoiseT ObserveNoiseT;                                                                                        \
                                                                                                                               \
  static_assert(is_vector<StateT>::value, "StateT must be an Eigen3 vector");                                                  \
  static_assert(is_vector<ObserveT>::value, "ObserveT must be an Eigen3 vector");                                              \
  static_assert(is_vector<ObserveNoiseT>::value, "ObserveNoiseT must be an Eigen3 vector");                                    \
                                                                                                                               \
  static_assert(is_same_scalar_type<StateT, ObserveT>::value, "StateT and ObserveT must have the same scalar type");           \
  static_assert(is_same_scalar_type<StateT, ObserveNoiseT>::value, "StateT and ObserveNoiseT must have the same scalar type"); \
                                                                                                                               \
  typedef typename MatrixTraits<StateT>::scalar_type scalar_type;                                                              \
                                                                                                                               \
  const static int _state_dim = MatrixTraits<StateT>::rows;                                                                    \
  const static int _observe_dim = MatrixTraits<ObserveT>::rows;                                                                \
  const static int _observe_noise_dim = MatrixTraits<ObserveNoiseT>::rows;                                                     \
                                                                                                                               \
  typedef Matrix<scalar_type, _observe_dim, _state_dim> ObserveStateGradT;                                                     \
  typedef Matrix<scalar_type, _observe_dim, _observe_noise_dim> ObserveNoiseGradT;
