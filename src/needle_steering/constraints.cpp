#include "needle_steering.hpp"

namespace Needle {

  PositionError::PositionError(LocalConfigurationPtr cfg, const Vector6d& target_pos, const Vector3d position_error_relax, NeedleProblemHelperPtr helper) : cfg(cfg), target_pose(expUp(target_pos)), position_error_relax(position_error_relax), body(cfg->GetBodies()[0]), helper(helper) {}

  VectorXd PositionError::operator()(const VectorXd& a) const {
    assert(a.size() == 6);
    //Matrix4d X = cfg->pose * expUp(a);
    //Vector6d x;
    //x.head<3>() = X.block<3, 1>(0, 3) - target_pos.head<3>();
    //x.tail<3>() = Vector3d::Zero();
    Matrix4d current_pose = cfg->pose * expUp(a);
    Vector3d orientation_error = logDown(current_pose.inverse() * target_pose).tail<3>();
    Vector3d position_error = (current_pose.block<3, 1>(0, 3) - target_pose.block<3, 1>(0, 3)).array().abs();
    position_error = (position_error - this->position_error_relax).cwiseMax(Vector3d::Zero());
    return concat(position_error, orientation_error);
  }

  PoseError::PoseError(LocalConfigurationPtr cfg0, LocalConfigurationPtr cfg1, NeedleProblemHelperPtr helper) : cfg0(cfg0), cfg1(cfg1), helper(helper) {}

  VectorXd PoseError::operator()(const VectorXd& a) const {
    Matrix4d pose1 = cfg0->pose * expUp(a.topRows(6));
    Matrix4d pose2 = cfg1->pose * expUp(a.middleRows(6,6));
    double Delta = a(12);
    double curvature_or_radius;
    switch (helper->curvature_constraint) {
      case NeedleProblemHelper::ConstantRadius:
        curvature_or_radius = 1.0 / helper->r_min;
        break;
      case NeedleProblemHelper::BoundedRadius:
        curvature_or_radius = a(13);
        break;
      SWITCH_DEFAULT;
    }
    double theta;
    switch (helper->curvature_formulation) {
      case NeedleProblemHelper::UseCurvature:
        theta = Delta * curvature_or_radius;
        break;
      case NeedleProblemHelper::UseRadius:
        theta = Delta / curvature_or_radius;
        break;
      SWITCH_DEFAULT;
    }
    Vector6d v; v << 0, 0, Delta, theta, 0, 0;
    return logDown((pose1 * expUp(v)).inverse() * pose2);
  }


  SpeedDeviationError::SpeedDeviationError(double deviation, NeedleProblemHelperPtr helper) : deviation(deviation), helper(helper) {}

  VectorXd SpeedDeviationError::operator()(const VectorXd& a) const {
    Vector1d x;
    x[0] = sqrt((a.array() - deviation).square().sum()) - deviation * 0.5;
    return x;
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
      case NeedleProblemHelper::Form2: {
        return logDown(helper->TransformPose(pose1, phi, Delta, curvature_or_radius).inverse() * pose2);
      }
      SWITCH_DEFAULT;
    }
  }

  int ControlError::outputSize() const {
    switch (helper->formulation) {
      case NeedleProblemHelper::Form1:
      case NeedleProblemHelper::Form2:
        return 6;
      SWITCH_DEFAULT;
    }
  }

}
