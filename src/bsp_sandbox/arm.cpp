#include "arm.hpp"

using namespace BSP;
using namespace Geometry2D;

namespace ArmBSP {

  Beam2D camera_fov(double x, double y, double angle, double camera_depth, double camera_span_angle) {
    Beam2D out;
    double x0 = x, y0 = y;
    double skew_depth = camera_depth / cos(camera_span_angle * 0.5);
    double x1 = x0 + skew_depth * cos(angle - camera_span_angle * 0.5),
        y1 = y0 + skew_depth * sin(angle - camera_span_angle * 0.5);
    double x2 = x0 + skew_depth * cos(angle + camera_span_angle * 0.5),
        y2 = y0 + skew_depth * sin(angle + camera_span_angle * 0.5);
    out.base = Vector2d(x0, y0);
    out.a = Vector2d(x1, y1);
    out.b = Vector2d(x2, y2);
    return out;
  }

  ArmBSPPlanner::ArmBSPPlanner() : BSPPlanner<ArmBSPProblemHelper>() {}

  void ArmBSPPlanner::initialize() {
    BSPPlanner<ArmBSPProblemHelper>::initialize(); 
    helper->link_lengths = link_lengths;
    link_total_lengths = helper->link_total_lengths = link_lengths.array().sum();
    helper->base_config = base_config;
    helper->real_object_pos = real_object_pos;
    helper->camera_base = camera_base;
    helper->camera_depth = camera_depth;
    helper->camera_span_angle = camera_span_angle;
    helper->n_fov_parts = n_fov_parts;
    //partition(camera, n_fov_parts, &helper->partitioned_fov);
  }

  Vector2d ArmBSPPlanner::get_feasible_pos(const Vector2d& pos) const {
    Vector2d base = base_config.head<2>();
    double distance = (pos - base).norm();
    if (distance > 0.8 * link_total_lengths) {
      return base + (pos - base) / distance * 0.8 * link_total_lengths;
    } else {
      return pos;
    }
  }

  void ArmBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt, bool is_first_time) {
    opt.max_iter_                   = 350;
    opt.merit_error_coeff_          = 100;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 2;
    opt.trust_shrink_ratio_         = 0.8;
    opt.trust_expand_ratio_         = 1.2;
    opt.min_trust_box_size_         = 0.001;
    opt.min_approx_improve_         = 1e-2;
    opt.min_approx_improve_frac_    = 1e-4;
    opt.improve_ratio_threshold_    = 0.25;
    opt.trust_box_size_             = 1;
    opt.cnt_tolerance_              = 1e-06;
  }

  void ArmBSPPlanner::custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state) {
    if (truncated_gaussian == WithTruncatedGaussian) {
      assert (state != NULL);
      assert (sigma != NULL);

      vector<Beam2D> actual_fov;
      helper->fov_from_state(actual_state, &actual_fov);

      if (!inside(real_object_pos, actual_fov)) { // truncate current belief if object is not in view
        vector<Beam2D> cur_fov;
        helper->fov_from_state(*state, &cur_fov);
        Vector2d new_state;
        Matrix2d new_sigma;
        truncate_belief(cur_fov, state->tail<2>(), sigma->bottomRightCorner<2, 2>(), &new_state, &new_sigma);
        state->tail<2>() = get_feasible_pos(new_state);
        sigma->bottomRightCorner<2, 2>() = ensure_precision(new_sigma);
      }
    }
  }
  
  ArmBSPProblemHelper::ArmBSPProblemHelper() : BSPProblemHelper<ArmBeliefFunc>() {
    input_dt = 1.0;
    set_state_dim(6);
    set_sigma_dof(21);
    set_observe_dim(6);
    set_control_dim(4);

    double state_lbs_array[] = {-INFINITY, -INFINITY, -INFINITY, PI/18, -15, -15};
    double state_ubs_array[] = {INFINITY, INFINITY, INFINITY, PI*17/18, 15, 15};
    double control_lbs_array[] = {-PI/32, -PI/32, -PI/32, -PI/32};
    double control_ubs_array[] = {PI/32, PI/32, PI/32, PI/32};

    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));

    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim) * 0.1);
  }

  void ArmBSPProblemHelper::initialize() {
    BSPProblemHelper<ArmBeliefFunc>::initialize();
    belief_constraints.clear();
  }

  void ArmBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    VectorOfVectorPtr f(new ArmGoalError(boost::static_pointer_cast<ArmBSPProblemHelper>(this->shared_from_this())));
    Vector2d coeffs = Vector2d::Ones();
    prob.addConstraint(ConstraintPtr(new ConstraintFromFunc(f, state_vars.row(T), coeffs, EQ, "goal")));
  }

  TransT ArmBSPProblemHelper::angle_to_transform(const Vector3d& angle) const {
    TransT mat(2, 3);
    double t0 = base_config(2) + angle(0);
    double t1 = t0 + angle(1);
    double t2 = t1 + angle(2);
    mat << cos(t0), cos(t1), cos(t2),
           sin(t0), sin(t1), sin(t2);
    return mat;
  }

  Vector8d ArmBSPProblemHelper::angle_to_endpoints(const Vector3d& angle) const {
    Vector3d l1; l1 << link_lengths(0), 0, 0;
    Vector3d l2; l2 << link_lengths(0), link_lengths(1), 0;
    Vector3d l3; l3 << link_lengths(0), link_lengths(1), link_lengths(2);
    TransT mat = angle_to_transform(angle);
    Vector8d res;
    res.segment<2>(0) = base_config.head<2>();
    res.segment<2>(2) = base_config.head<2>() + mat * l1;
    res.segment<2>(4) = base_config.head<2>() + mat * l2;
    res.segment<2>(6) = base_config.head<2>() + mat * l3;
    return res;
  }

  vector<Segment> ArmBSPProblemHelper::angle_to_segments(const Vector3d& angle) const {
    vector<Segment> segs;
    Vector8d endpoints = angle_to_endpoints(angle);
    Segment seg1; seg1.a = endpoints.segment<2>(0); seg1.b = endpoints.segment<2>(2);
    Segment seg2; seg2.a = endpoints.segment<2>(2); seg2.b = endpoints.segment<2>(4);
    Segment seg3; seg3.a = endpoints.segment<2>(4); seg3.b = endpoints.segment<2>(6);
    segs.push_back(seg1); segs.push_back(seg2); segs.push_back(seg3);
    return segs;
  }

  Vector2d ArmBSPProblemHelper::angle_to_pos(const Vector3d& angle) const {
    Vector3d l3; l3 << link_lengths(0), link_lengths(1), link_lengths(2);
    TransT mat = angle_to_transform(angle);
    return mat * l3 + base_config.head<2>();
  }

  void ArmBSPProblemHelper::fov_from_state(const StateT& x, vector<Beam2D>* out_fov) const {
    Beam2D fov = camera_fov(camera_base.x(), camera_base.y(), x(3), camera_depth, camera_span_angle);
    vector<Beam2D> partitioned_fov;
    partition(fov, n_fov_parts, &partitioned_fov);
    truncate_beams(partitioned_fov, angle_to_segments(x.head<3>()), out_fov);
  }

  ArmStateFunc::ArmStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ArmStateFunc::ArmStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  StateT ArmStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT ret;
    ret.head<4>() = x.head<4>() + u + 0.01 * m.head<4>();
    ret.tail<2>() = x.tail<2>() + 0.0001 * m.tail<2>();
    return ret;
  }

  ArmObserveFunc::ArmObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ArmObserveFunc::ArmObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  ObserveT ArmObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return x + 0.01 * n;
  }

  ObserveT ArmObserveFunc::observe_masks_from_object_position(const StateT& x, const Vector2d& object_pos, double approx_factor) const {
    vector<Beam2D> cur_fov;
    arm_helper->fov_from_state(x, &cur_fov);

    ObserveT ret(observe_dim);
    ret(0) = ret(1) = ret(2) = ret(3) = 1;
    if (approx_factor < 0) {
      if (inside(object_pos, cur_fov)) {
        ret(4) = ret(5) = 1;
      } else {
        ret(4) = ret(5) = 0;
      }
    } else {
      double dist = Geometry2D::sgndist(x.tail<2>(), cur_fov);
      double tol = 0.1;
      ret(4) = ret(5) = 1 - sigmoid(approx_factor * (dist + tol));
    }
    return ret;
  }

  ObserveT ArmObserveFunc::observation_masks(const StateT& x, double approx_factor) const {
    return observe_masks_from_object_position(x, x.tail<2>(), approx_factor);
  }

  ObserveT ArmObserveFunc::real_observation(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT ret;
    ret.head<4>() = x.head<4>() + 0.01 * n.head<4>();
    ret.tail<2>() = arm_helper->real_object_pos + 0.01 * n.tail<2>();
    return ret;
  }

  ObserveT ArmObserveFunc::real_observation_masks(const StateT& x) const {
    return observe_masks_from_object_position(x, arm_helper->real_object_pos, -1);
  }

  ArmBeliefFunc::ArmBeliefFunc() : BeliefFunc<ArmStateFunc, ArmObserveFunc, BeliefT>() {
    set_approx_factor(0.5);
  }

  ArmBeliefFunc::ArmBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
    BeliefFunc<ArmStateFunc, ArmObserveFunc, BeliefT>(helper, f, h), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {
    set_approx_factor(0.5);
  }

  ArmGoalError::ArmGoalError(ArmBSPProblemHelperPtr helper) : helper(helper) {}

  VectorXd ArmGoalError::operator()(const VectorXd& a) const {
    assert (a.size() == 5);
    return a.tail<2>() - helper->angle_to_pos(a.head<3>());
  }

  ArmPlotter::ArmPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  void ArmPlotter::draw_robot(const Vector3d& x, QPainter& painter) {
    Vector8d endpoints = arm_helper->angle_to_endpoints(x);
    draw_line(endpoints(0), endpoints(1), endpoints(2), endpoints(3), painter);
    draw_line(endpoints(2), endpoints(3), endpoints(4), endpoints(5), painter);
    draw_line(endpoints(4), endpoints(5), endpoints(6), endpoints(7), painter);
  }

  void ArmPlotter::draw_beam_2d(const Beam2D& beam, QPainter& painter) {
    QPolygonF polygon;
    polygon << QPointF(scale_x(beam.base.x()), scale_y(beam.base.y()));
    polygon << QPointF(scale_x(beam.a.x()), scale_y(beam.a.y()));
    polygon << QPointF(scale_x(beam.b.x()), scale_y(beam.b.y()));
    painter.drawPolygon(polygon);
  }

  void ArmPlotter::draw_beams(const vector<Beam2D>& beams, QPainter& painter) {
    if (beams.size() <= 0) {
      return;
    }
    QPolygonF polygon;
    polygon << QPointF(scale_x(beams[0].base.x()), scale_y(beams[0].base.y()));
    for (int i = 0; i < beams.size(); ++i) {
      polygon << QPointF(scale_x(beams[i].a.x()), scale_y(beams[i].a.y()));
      polygon << QPointF(scale_x(beams[i].b.x()), scale_y(beams[i].b.y()));
    }
    painter.drawPolygon(polygon);
  }

  void ArmPlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
    assert(distmap != NULL);
    if (state != NULL) {
      vector<Beam2D> cur_fov;
      arm_helper->fov_from_state(*state, &cur_fov);
      for (int j = 0; j < distmap->height(); ++j) {
        QRgb *line = (QRgb*) distmap->scanLine(j);
        for (int i = 0; i < distmap->width(); ++i) {
          double x = unscale_x(i),
                 y = unscale_y(j);
          double dist = sgndist(Vector2d(x, y), cur_fov);
          double grayscale;
          if (approx_factor > 0) {
            grayscale = 1./(1. + exp(approx_factor*dist));
          } else {
            grayscale = dist <= 0 ? 1 : 0;
          }
          line[i] = qRgb(grayscale*255, grayscale*255, grayscale*255);
        }
      }
    }
  }

  ArmOptPlotter::ArmOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   ArmPlotter(x_min, x_max, y_min, y_max, helper, parent),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void ArmOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      if (states.size() > 0) {
        compute_distmap(&distmap, &states[0], arm_helper->belief_func->approx_factor);
      } else {
        compute_distmap(&distmap, NULL, arm_helper->belief_func->approx_factor);
      }
    }
    if (states.size() > 0) {
      //painter.drawImage(0, 0, distmap);
      QPen robot_start_pen(QColor(255, 0, 0, 255), 8, Qt::SolidLine);
      QPen robot_middle_pen(QColor(255, 0, 0, 100), 8, Qt::SolidLine);
      QPen robot_end_pen(QColor(255, 0, 0, 255), 8, Qt::SolidLine);

      painter.setRenderHint(QPainter::Antialiasing);
      painter.setRenderHint(QPainter::HighQualityAntialiasing);
      painter.setPen(robot_start_pen);
      draw_robot(states[0].head<3>(), painter);

      painter.setPen(robot_middle_pen);
      for (int i = 1; i < states.size() - 1; ++i) {
        draw_robot(states[i].head<3>(), painter);
      }

      if (states.size() > 1) {
        painter.setPen(robot_end_pen);
        draw_robot(states[states.size() - 1].head<3>(), painter);
      }

      painter.setPen(QPen(QColor(255, 255, 255, 10), 0, Qt::SolidLine));
      {
        QBrush prev_brush = painter.brush();
        painter.setBrush(QBrush(QColor(255, 255, 255, 20)));
        for (int i = 0; i < states.size(); ++i) {
          //cout << "angle: " << states[i](3) << endl;
          vector<Beam2D> fov;
          arm_helper->fov_from_state(states[i], &fov);
          //Beam2D fov = camera_fov(arm_helper->camera_base.x(), arm_helper->camera_base.y(), states[i](3), arm_helper->camera_depth, arm_helper->camera_span_angle);
          //for (int j = 0; j < fov.size(); ++j) {
            draw_beams(fov, painter);
          //}
        }
        painter.setBrush(prev_brush);
      }
      if (states.size() > 0) {
        {
          QBrush prev_brush = painter.brush();
          painter.setBrush(QBrush(QColor(255, 215, 0)));
          painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));
          draw_point_with_border(states[0](4), states[0](5), 0.2, 0.2, painter);
        }
        {
          QBrush prev_brush = painter.brush();
          painter.setBrush(QBrush(Qt::green));
          painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));
          draw_point_with_border(arm_helper->real_object_pos.x(), arm_helper->real_object_pos.y(), 0.3, 0.3, painter);
          painter.setBrush(prev_brush);
        }
      }
    }
  }

  void ArmOptPlotter::update_plot_data(void* data) {
    DblVec* xvec = (DblVec* ) data;
    vector<StateT> new_states;
    vector<VarianceT> new_sigmas;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = arm_helper->belief_func->approx_factor;
    BeliefT cur_belief;
    arm_helper->belief_func->compose_belief(arm_helper->start, matrix_sqrt(arm_helper->start_sigma), &cur_belief);
    for (int i = 0; i <= arm_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      arm_helper->belief_func->extract_state(cur_belief, &cur_state);
      arm_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < arm_helper->T) cur_belief = arm_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, arm_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;
    this->repaint();
  }

  ArmSimulationPlotter::ArmSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   ArmPlotter(x_min, x_max, y_min, y_max, helper, parent),
   distmap(400, 400, QImage::Format_RGB32) {}

  void ArmSimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    distmap = QImage(width(), height(), QImage::Format_RGB32);
    if (simulated_positions.size() > 0) {
      compute_distmap(&distmap, &simulated_positions.back(), -1);
    } else {
      compute_distmap(&distmap, NULL, -1);
    }
    painter.drawImage(0, 0, distmap);
    if (simulated_positions.size() <= 0) {
      return;
    }
    QPen sim_prev_robot_pen(QColor(255, 0, 0, 100), 8, Qt::SolidLine);
    QPen sim_cur_robot_pen(QColor(255, 0, 0, 255), 8, Qt::SolidLine);
    QPen real_object_pos_pen(Qt::green, 8, Qt::SolidLine);
    QPen object_position_belief_pen(QColor(255,215,0), 4, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(sim_prev_robot_pen);
    for (int i = 0; i < (int) simulated_positions.size() - 1; ++i) {
      draw_robot(simulated_positions[i].head<3>(), painter);
    }
    painter.setPen(sim_cur_robot_pen);
    if (simulated_positions.size() > 0) {
      draw_robot(simulated_positions.back().head<3>(), painter);
    }
    {
      QBrush prev_brush = painter.brush();
      painter.setBrush(QBrush(Qt::green));
      painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));
      draw_point_with_border(arm_helper->real_object_pos.x(), arm_helper->real_object_pos.y(), 0.3, 0.3, painter);
      painter.setBrush(prev_brush);
    }
    //painter.setPen(real_object_pos_pen);
    {
      QBrush prev_brush = painter.brush();
      painter.setBrush(QBrush(QColor(255, 215, 0)));
      painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));

      //painter.setPen(object_position_belief_pen);
      for (int i = 0; i < object_position_beliefs.size(); ++i) {
        draw_point_with_border(object_position_beliefs[i].x(), object_position_beliefs[i].y(), 0.2, 0.2, painter);
      }
      painter.setBrush(prev_brush);
    }
    painter.setPen(object_position_belief_pen);
      //for (int i = 0; i < object_position_sigmas.size(); ++i) {
        draw_ellipse(object_position_beliefs.back(), object_position_sigmas.back(), painter);
      //  cout << "position: " << object_position_beliefs[i].transpose() << endl;
      //  cout << "sigma: " << object_position_sigmas[i] << endl;
      //}
  }

  void ArmSimulationPlotter::update_plot_data(void* data_x, void* data_sim) {
    simulated_positions = *((vector<StateT>*) data_sim);
    object_position_beliefs.push_back(arm_helper->start.tail<2>());
    object_position_sigmas.push_back(arm_helper->start_sigma.bottomRightCorner<2, 2>());
    //DblVec* xvec = (DblVec* ) data_x;
    //if (xvec->size() > 0) {
    //  //vector<StateT> new_states;
    //  //vector<VarianceT> new_sigmas;
    //  //BeliefT cur_belief;
    //  //arm_helper->belief_func->compose_belief(arm_helper->start, matrix_sqrt(arm_helper->start_sigma), &cur_belief);
    //  //for (int i = 0; i <= arm_helper->T; ++i) {
    //  //  StateT cur_state;
    //  //  VarianceT cur_sigma;
    //  //  arm_helper->belief_func->extract_state(cur_belief, &cur_state);
    //  //  arm_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
    //  //  new_states.push_back(cur_state);
    //  //  new_sigmas.push_back(cur_sigma);
    //  //  if (i < arm_helper->T) cur_belief = arm_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, arm_helper->control_vars.row(i)));
    //  //}
    //  //states = new_states;
    //  //sigmas = new_sigmas;

    //  
    //} else {
    //  states.clear();
    //  states.push_back(arm_helper->start);
    //  sigmas.clear();
    //  sigmas.push_back(arm_helper->start_sigma);
    //  //object_position_beliefs.push_back(new_states[0].tail<2>());
    //  //object_position_sigmas.push_back(new_sigmas[0].bottomRightCorner<2, 2>());
    //}

    this->repaint();
  }

  ArmOptimizerTask::ArmOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  ArmOptimizerTask::ArmOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void ArmOptimizerTask::stage_plot_callback(boost::shared_ptr<ArmOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
    //wait_to_proceed(boost::bind(&ArmOptPlotter::update_plot_data, plotter, &x));
  }

  void ArmOptimizerTask::run() {
    srand(static_cast<unsigned int>(std::time(0)));
    int T = 12;
    bool plotting = false;
    bool first_step_only = false;
    int method = 2;
    int truncated_gaussian = 1;
    double noise_level = 0.2;
    double base_vec_array[] = {0, -1, 0};
    double start_vec_array[] = {PI/4+PI/16, PI/2, PI/4+PI/16, PI/2, 5, 5};
    vector<double> base_vec(base_vec_array, end(base_vec_array));
    vector<double> start_vec(start_vec_array, end(start_vec_array));
    {
      Config config;
      config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter<int>("method", &method, "method"));
      config.add(new Parameter<int>("truncated_gaussian", &truncated_gaussian, "truncated_gaussian"));
      config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
      config.add(new Parameter<bool>("first_step_only", &first_step_only, "first_step_only"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }
    Vector6d start = toVectorXd(start_vec);
    Vector3d base_config = toVectorXd(base_vec);
    Matrix6d start_sigma = Matrix6d::Identity() * 0.01;
    start_sigma.bottomRightCorner<2, 2>() = Matrix2d::Identity() * 20;

    deque<Vector4d> initial_controls;
    for (int i = 0; i < T; ++i) {
      initial_controls.push_back(Vector4d::Zero());
    }

    Vector3d link_lengths; link_lengths << 5, 5, 4;

    //Beam2D camera;
    //camera.base = Vector2d(0, 0);
    //camera.a = Point(5, 10);
    //camera.b = Point(-5, 10);

    ArmBSPPlannerPtr planner(new ArmBSPPlanner());

    planner->start = start;
    planner->base_config = base_config;
    planner->start_sigma = start_sigma;
    planner->link_lengths = link_lengths;
    planner->T = T;
    planner->noise_level = noise_level;
    planner->method = method;
    planner->truncated_gaussian = truncated_gaussian;
    planner->controls = initial_controls;
    planner->camera_base = Vector2d(0, 0);
    planner->camera_span_angle = PI / 4;
    planner->camera_depth = 10;
    planner->n_fov_parts = 30;
    planner->real_object_pos = Vector2d(6, 6);
    planner->initialize();

    boost::shared_ptr<ArmSimulationPlotter> sim_plotter;
    boost::shared_ptr<ArmOptPlotter> opt_plotter;
    if (plotting) {
      double x_min = -12.4, x_max = 12.4, y_min = -2, y_max = 12;
      sim_plotter.reset(create_plotter<ArmSimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      opt_plotter.reset(create_plotter<ArmOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      opt_plotter->show();
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    if (plotting) {
      //opt_callback = boost::bind(&ArmOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    if (plotting) {
      emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions, false);
    }

    Vector6d state_error = Vector6d::Ones() * 10000;

    while (!planner->finished()) {
      if (state_error.array().abs().sum() < 0.0001) {
        // nothing
      } else if (state_error.array().abs().sum() < 0.10) {
        planner->solve(opt_callback, 9, 1);
      } else {
        planner->solve(opt_callback, 1, 3);
      }
      if (first_step_only) break;
      state_error = planner->simulate_executions(1);
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions, false);
        //sim_plotter->update_plot_data(&planner->result, &planner->simulated_positions);
      }
    }

    emit finished_signal();
  }
}

using namespace ArmBSP;

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  ArmOptimizerTask* task = new ArmOptimizerTask(argc, argv, &app);
  QTimer::singleShot(0, task, SLOT(run_slot()));
  QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
  return app.exec();
}
