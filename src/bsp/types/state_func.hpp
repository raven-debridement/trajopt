#pragma once

#define define_state_func_types()                                                                                                  \
  typedef _StateT StateT;                                                                                                          \
  typedef _ControlT ControlT;                                                                                                      \
  typedef _StateNoiseT StateNoiseT;                                                                                                \
                                                                                                                                   \
  static_assert(is_vector<StateT>::value, "StateT must be an Eigen3 vector");                                                      \
  static_assert(is_vector<ControlT>::value, "ControlT must be an Eigen3 vector");                                                  \
  static_assert(is_vector<StateNoiseT>::value, "StateNoiseT must be an Eigen3 vector");                                            \
                                                                                                                                   \
  static_assert(is_same_scalar_type<StateT, ControlT>::value, "StateT and ControlT must have the same scalar type");               \
  static_assert(is_same_scalar_type<StateT, StateNoiseT>::value, "StateT and StateNoiseT vectors must have the same scalar type"); \
                                                                                                                                   \
  typedef typename MatrixTraits<StateT>::scalar_type scalar_type;                                                                  \
                                                                                                                                   \
  const static int _state_dim = MatrixTraits<StateT>::rows;                                                                        \
  const static int _control_dim = MatrixTraits<ControlT>::rows;                                                                    \
  const static int _state_noise_dim = MatrixTraits<StateNoiseT>::rows;                                                             \
                                                                                                                                   \
  typedef Matrix<scalar_type, _state_dim, _state_dim> StateGradT;                                                                  \
  typedef Matrix<scalar_type, _state_dim, _control_dim> ControlGradT;                                                              \
  typedef Matrix<scalar_type, _state_dim, _state_noise_dim> StateNoiseGradT;
