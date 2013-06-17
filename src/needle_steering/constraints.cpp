#include "needle_steering.hpp"

namespace Needle {

  PositionError::PositionError(LocalConfigurationPtr cfg, const VectorXd& target_pos) : cfg(cfg), target_pos(target_pos), body(cfg->GetBodies()[0]) {}

  VectorXd PositionError::operator()(const VectorXd& a) const {
    Matrix4d X = cfg->pose * expUp(a);
    VectorXd x(6);
    x.head<3>() = X.block<3, 1>(0, 3) - target_pos.head<3>();
    x.tail<3>() = Vector3d::Zero();
    return logDown((cfg->pose * expUp(a)).inverse() * expUp(target_pos));
  }

  ControlError::ControlError(LocalConfigurationPtr cfg0, LocalConfigurationPtr cfg1, NeedleProblemHelperPtr helper) : cfg0(cfg0), cfg1(cfg1), body(cfg0->GetBodies()[0]), helper(helper) {}

  VectorXd ControlError::operator()(const VectorXd& a) const {
    Matrix4d pose1 = cfg0->pose * expUp(a.topRows(6));
    Matrix4d pose2 = cfg1->pose * expUp(a.middleRows(6,6));
    double phi = a(12), Delta = a(13);
    double curvature_or_radius;
    switch (helper->curvature_formulation) {
      case NeedleProblemHelper::UseCurvature:
        switch (helper->curvature_constraint) {
          case NeedleProblemHelper::ConstantRadius:
            curvature_or_radius = 1.0 / helper->r_min;
            break;
          case NeedleProblemHelper::BoundedRadius:
            curvature_or_radius = a(14);
            break;
          SWITCH_DEFAULT;
        }
        break;
      case NeedleProblemHelper::UseRadius:
        switch (helper->curvature_constraint) {
          case NeedleProblemHelper::ConstantRadius:
            curvature_or_radius = helper->r_min;
            break;
          case NeedleProblemHelper::BoundedRadius:
            curvature_or_radius = a(14);
            break;
          SWITCH_DEFAULT;
        }
        break;
      SWITCH_DEFAULT;
    }
    switch (helper->formulation) {
      case NeedleProblemHelper::Form1:
      case NeedleProblemHelper::Form3: {
        return logDown(helper->TransformPose(pose1, phi, Delta, curvature_or_radius).inverse() * pose2);
      }
      case NeedleProblemHelper::Form2: {
        VectorXd trans1 = pose1.block<3, 1>(0, 3);
        VectorXd trans2 = pose2.block<3, 1>(0, 3);
        Matrix3d R1 = pose1.block<3, 3>(0, 0);
        Vector3d xyz = R1.transpose() * (trans2 - trans1);
        double x = xyz[0], y = xyz[1], z = xyz[2];
        Vector3d err;// err <<
        //  phi - atan2(x, -y),
        //  bound_inf(radius - (x*x + y*y + z*z) / (2 * sqrt(x*x + y*y)), 100), // sketchy lol
        //  Delta - radius * atan2(z, radius - sqrt(x*x + y*y));
        return err;
      }
      SWITCH_DEFAULT;
    }
  }

  int ControlError::outputSize() const {
    switch (helper->formulation) {
      case NeedleProblemHelper::Form1:
        return 6;
      case NeedleProblemHelper::Form2:
        return 3;
      case NeedleProblemHelper::Form3:
        return 6;
      SWITCH_DEFAULT;
    }
  }

}
