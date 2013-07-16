#pragma once
#include "bsp/bsp.hpp"
#include <boost/math/distributions.hpp>

namespace Geometry2D {
  const double eps = 1e-5;
  const double m_sqrt_pi_2 = sqrt(PI / 2.0);
  const double m_sqrt_1_2 = sqrt(0.5);
  const boost::math::normal_distribution<> standard_normal;
  typedef Vector2d Point;

  // Assume that the points {base, a, b} are ordered counterclockwise
  struct Beam2D {
    Point base;
    Point a;
    Point b;
    friend ostream& operator<< (ostream &out, const Beam2D& beam);
  };

  struct Segment {
    Point a;
    Point b;
  };

  double dist_to_segment(const Point& p, const Point& v, const Point& w);

  bool line_intersection(const Point& p1, const Point& p2, const Point& p3, const Point& p4, Point* intersection);

  bool segment_intersects_segment(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

  Beam2D truncate_beam(const Beam2D& beam, const Segment& segment);

  Beam2D truncate_beam(const Beam2D& beam, const vector<Segment>& segments);

  void truncate_beams(const vector<Beam2D> beams, const vector<Segment>& segments, vector<Beam2D>* truncated_beams);

  bool inside(const Point& p, const Beam2D& beam);

  bool inside(const Point& p, const vector<Beam2D>& beams);

  double sgndist(const Point& p, const Beam2D& beam);

  double sgndist(const Point& p, const vector<Beam2D>& beams);

  void partition(const Beam2D& beam, int nparts, vector<Beam2D>* beams);

  template< class VecT >
  Matrix< typename BSP::MatrixTraits<VecT>::scalar_type, BSP::MatrixTraits<VecT>::rows, BSP::MatrixTraits<VecT>::rows> // return matrix type
  rotation_mat(const VecT& a, const VecT& v) {
    typedef Matrix< typename BSP::MatrixTraits<VecT>::scalar_type, BSP::MatrixTraits<VecT>::rows, BSP::MatrixTraits<VecT>::rows> MatT;
    assert (a.size() == v.size());
    if ((a - v).norm() < eps) {
      return MatT::Identity(a.size(), a.size());
    }
    VecT b = v - a * a.dot(v);
    double norm_b = b.norm();
    if (norm_b < eps) {
      return -MatT::Identity(a.size(), a.size());
    }
    b /= norm_b;

    double v_a = v.dot(a);
    double v_b = v.dot(b);

    MatT R(a.size(), a.size());

    for (int i = 0; i < a.size(); ++i) {
      VecT x = VecT::Zero(a.size());
      x(i) = 1.0;
      double x_a = x.dot(a);
      double x_b = x.dot(b);
      VecT x_res = x - x_a*a - x_b*b;
      x = x_res + (v_a*x_a - v_b*x_b)*a + (v_b*x_a + v_a*x_b)*b;
      R.col(i) = x;
    }
    return R;
  }

  void truncate_univariate_gaussian(double x, double cur_mean, double cur_var, double *out_mean, double *out_var);

  //void truncate_standard_gaussian(double *mean, double *var, double z);
  
  // truncate to cx <= d
  template< class VecT, class MatT >
  void truncate_gaussian(const VecT& c, double d, const VecT& cur_mean, const MatT& cur_cov, VecT* out_mean, MatT* out_cov) {
    assert (out_mean != nullptr);
    assert (out_cov != nullptr);
    assert (cur_mean.size() == cur_cov.rows());
    assert (cur_mean.size() == cur_cov.cols());
    double y_mean = c.dot(cur_mean);
    double y_var = c.transpose() * cur_cov * c;
    double y_new_mean, y_new_var;
    truncate_univariate_gaussian(d, y_mean, y_var, &y_new_mean, &y_new_var);
    VecT xy_cov = cur_cov * c;
    VecT L = xy_cov / y_var;
    *out_mean = cur_mean + L * (y_new_mean - y_mean);
    *out_cov = cur_cov + (y_new_var / y_var - 1.0) * (L * xy_cov.transpose());
  }

  //template< class VecT, class MatT >
  //void truncate_gaussian(const VecT& c, double d, const VecT& cur_mean, const MatT& cur_cov, VecT* out_mean, MatT* out_cov) {
  //  assert (out_mean != nullptr);
  //  assert (out_cov != nullptr);
  //  assert (cur_mean.size() == cur_cov.rows());
  //  assert (cur_mean.size() == cur_cov.cols());

  //  LLT<MatT> llt(cur_cov);
  //  MatT M = llt.matrixU();
  //  VecT v = M*c;
  //  double norm_v = v.norm();
  //  if (norm_v < eps) {
  //    cout << "warning: degenerate truncation (something is probably wrong)" << endl;
  //    *out_mean = cur_mean;
  //    *out_cov = cur_cov;
  //    return;
  //  }
  //  double z = (c.dot(cur_mean) - d) / norm_v;
  //  v /= norm_v;
  //  VecT e = VecT::Zero(cur_mean.size()); e(0) = 1.0;
  //  MatT R = rotation_mat(e, v);
  //  double mean, var;
  //  truncate_standard_gaussian(&mean, &var, -z);
  //  VecT b = VecT::Zero(cur_mean.size()); b(0) = mean;
  //  MatT B = MatT::Identity(cur_mean.size(), cur_mean.size()); B(0, 0) = var;
  //  *out_mean = M.transpose() * R * b + cur_mean;
  //  *out_cov = M.transpose() * R * B * R.transpose() * M;
  //}

  void truncate_belief(const vector<Beam2D>& beams, const Vector2d& cur_mean, const Matrix2d& cur_cov, Vector2d* out_mean, Matrix2d* out_cov);
}
