#include "needle_steering.hpp"

namespace Needle {

  inline double bound_inf(double result, double bound) {
    return min(max(result, -bound), bound);
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<double>& lbs, const vector<double>& ubs, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    int n_arr = name_prefix.size();
    assert(n_arr == newvars.size());

    vector<MatrixXi> index(n_arr);
    for (int i=0; i < n_arr; ++i) {
      newvars[i]->resize(rows, cols[i]);
      index[i].resize(rows, cols[i]);
    }

    vector<string> names;
    vector<double> all_lbs;
    vector<double> all_ubs;
    int var_idx = prob.getNumVars();
    for (int i=0; i < rows; ++i) {
      for (int k=0; k < n_arr; ++k) {
        for (int j=0; j < cols[k]; ++j) {
          index[k](i,j) = var_idx;
          names.push_back( (boost::format("%s_%i_%i")%name_prefix[k]%i%j).str() );
          all_lbs.push_back(lbs[k]);
          all_ubs.push_back(ubs[k]);
          ++var_idx;
        }
      }
    }
    prob.createVariables(names, all_lbs, all_ubs); // note that w,r, are both unbounded

    const vector<Var>& vars = prob.getVars();
    for (int k=0; k < n_arr; ++k) {
      for (int i=0; i < rows; ++i) {
        for (int j=0; j < cols[k]; ++j) {
          (*newvars[k])(i,j) = vars[index[k](i,j)];
        }
      }
    }
  }

  void AddVarArrays(OptProb& prob, int rows, const vector<int>& cols, const vector<string>& name_prefix, const vector<VarArray*>& newvars) {
    vector<double> lbs(newvars.size(), -INFINITY);
    vector<double> ubs(newvars.size(), INFINITY);
    AddVarArrays(prob, rows, cols, lbs, ubs, name_prefix, newvars);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, double lb, double ub, const string& name_prefix, VarArray& newvars) {
    vector<VarArray*> arrs(1, &newvars);
    vector<string> prefixes(1, name_prefix);
    vector<int> colss(1, cols);
    vector<double> lbs(1, lb);
    vector<double> ubs(1, ub);
    AddVarArrays(prob, rows, colss, lbs, ubs, prefixes, arrs);
  }

  void AddVarArray(OptProb& prob, int rows, int cols, const string& name_prefix, VarArray& newvars) {
    AddVarArray(prob, rows, cols, -INFINITY, INFINITY, name_prefix, newvars);
  }

  Matrix3d rotMat(const Vector3d& x) {
    Matrix3d out;
    out << 0, -x(2), x(1),
           x(2), 0, -x(0),
           -x(1), x(0), 0;
    return out;
  }

  Vector3d rotVec(const Matrix3d& X) {
    Vector3d out;
    out << X(2, 1), X(0, 2), X(1, 0);
    return out;
  }

  Matrix3d expA(const Vector3d& w) {
    double theta = w.norm();
    if (fabs(theta) < 1e-10) {
      return Matrix3d::Identity();
    } else {
      Matrix3d w_hat = rotMat(w);
      return Matrix3d::Identity() + w_hat / (theta*theta) * (1 - cos(theta)) + w_hat*w_hat / (theta*theta*theta) * (theta - sin(theta));
    }
  }

  Matrix3d logInvA(const Vector3d& w) {
    double theta = w.norm();
    Matrix3d w_hat = rotMat(w);
    if (fabs(theta) < 1e-8) {
      return Matrix3d::Identity();
    }
    return Matrix3d::Identity() - 0.5*w_hat + (2*sin(theta) - theta*(1 + cos(theta))) / (2 * theta*theta * sin(theta)) * w_hat*w_hat;
  }

  Matrix3d expRot(const Vector3d& x) {
    double rr = x.squaredNorm();
    if (fabs(rr) < 1e-10) {
      return Matrix3d::Identity();
    } else {
      double r = sqrt(rr);
      return rotMat(x * (sin(r) / r)) + Matrix3d::Identity() * cos(r) + (x*x.transpose()) * ((1 - cos(r)) / rr);
    }
  }

  Vector3d logRot(const Matrix3d& X) {
    // Using the old implementation since it seems more robust in practice
    Vector3d x;
    x << X(2, 1) - X(1, 2),
         X(0, 2) - X(2, 0),
         X(1, 0) - X(0, 1);
    double r = x.norm();
    double t = X(0, 0) + X(1, 1) + X(2, 2) - 1;

    if (fabs(r) < 1e-8) {
      return Vector3d::Zero();
    } else {
      return x * (atan2(r, t) / r);
    }
  }

  Matrix4d expUp(const Vector6d& x) {
    assert(x.size() == 6);
    Matrix4d X = Matrix4d::Identity();
    X.block<3, 3>(0, 0) = expRot(x.tail<3>());
    X.block<3, 1>(0, 3) = expA(x.tail<3>()) * x.head<3>();
    X(3, 3) = 1;
    return X;
  }

  Vector6d logDown(const Matrix4d& X) {
    VectorXd x(6);
    x.tail<3>() = logRot(X.block<3, 3>(0, 0));
    x.head<3>() = (expA(x.tail<3>())).inverse() * X.block<3, 1>(0, 3);
    return x;
  }

  OpenRAVE::Transform matrixToTransform(const Matrix4d& X) {
    OpenRAVE::TransformMatrix M;
    M.trans.x = X(0, 3);
    M.trans.y = X(1, 3);
    M.trans.z = X(2, 3);
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        M.m[row*4+col] = X(row, col);
      }
    }
    return OpenRAVE::Transform(M);
  }

  OpenRAVE::Transform vecToTransform(const Vector6d& x) {
    OpenRAVE::Transform T;
    OpenRAVE::Vector trans(x[0], x[1], x[2]);
    OpenRAVE::Vector rot(x[3], x[4], x[5]);
    T.trans = trans;
    T.rot = OpenRAVE::geometry::quatFromAxisAngle(rot);
    return T;
  }
}
