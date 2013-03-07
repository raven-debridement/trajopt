#pragma once

#include <Eigen/Core>
#include <vector>

namespace util {

template<typename T> inline std::vector<T> toDblVec(const Eigen::Matrix<T,Eigen::Dynamic,1>& x) {
  return std::vector<T>(x.data(), x.data()+x.size());
}
inline std::vector<double> toDblVec(const Eigen::Matrix<double,Eigen::Dynamic,1>& x) {
  return std::vector<double>(x.data(), x.data()+x.size());
}
template<typename T> inline Eigen::Matrix<T,Eigen::Dynamic,1> toVectorXd(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1> >(x.data(), x.size());
}
inline Eigen::VectorXd toVectorXd(const std::vector<double>& x) {
  return Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
}

template<typename T> inline BasicArray<T> toBasicArray(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& x) {
  return BasicArray<T>(x.rows(), x.cols(), x.data());
}
template<typename T> inline Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> toMatrixX(const BasicArray<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(x.data(), x.rows(), x.cols());
}

}

