#include "arm.hpp"

using namespace BSP;
using namespace Geometry2D;

namespace ArmBSP {

  ArmBSPPlanner::ArmBSPPlanner() : BSPPlanner<ArmBSPProblemHelper>() {}

  void ArmBSPPlanner::initialize() {
    BSPPlanner<ArmBSPProblemHelper>::initialize(); 
    helper->link_lengths = link_lengths;
    helper->base_config = base_config;
    helper->real_object_pos = real_object_pos;
    helper->camera = camera;
    helper->n_fov_parts = n_fov_parts;
    partition(camera, n_fov_parts, &helper->partitioned_fov);
  }

  Vector2d get_pos_within_region(const Vector2d& pos) {
    Vector2d out = pos;
    out.y() = fmax(fmin(out.y(), 10.), 0.);
    if (out.x() < 0) {
      out.x() = fmax(out.x(), -out.y() * 0.5);
    } else {
      out.x() = fmin(out.x(), out.y() * 0.5);
    }
    return out;
  }

  void ArmBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt) {
    opt.max_iter_                   = 250;
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
    assert (state != NULL);
    assert (sigma != NULL);
    vector<Beam2D> cur_fov;
    truncate_beams(helper->partitioned_fov, helper->angle_to_segments(state->head<3>()), &cur_fov);
    cout << "state: " << state->transpose() << endl;
    if (!inside(real_object_pos, cur_fov)) { // truncate current belief if object is not in view
      Vector2d new_state;
      Matrix2d new_sigma;
      truncate_belief(cur_fov, state->tail<2>(), sigma->bottomRightCorner<2, 2>(), &new_state, &new_sigma);
      state->tail<2>() = get_pos_within_region(new_state);
      sigma->bottomRightCorner<2, 2>() = ensure_precision(new_sigma);
    }
  }
  
  ArmBSPProblemHelper::ArmBSPProblemHelper() : BSPProblemHelper<ArmBeliefFunc>() {
    input_dt = 1.0;
    set_state_dim(5);
    set_sigma_dof(15);
    set_observe_dim(5);
    set_control_dim(3);
    set_state_bounds(concat(DblVec(3, -INFINITY), DblVec(2, -15)), concat(DblVec(3, INFINITY), DblVec(2, 15)));
    set_control_bounds(DblVec(3, -PI/32), DblVec(3, PI/32));
    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 100);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim) * 0.1);
  }

  void ArmBSPProblemHelper::initialize() {
    BSPProblemHelper<ArmBeliefFunc>::initialize();
    belief_constraints.clear();
    truncate_beams(partitioned_fov, angle_to_segments(start.head<3>()), &cur_fov);
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

  ArmStateFunc::ArmStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ArmStateFunc::ArmStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  StateT ArmStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT ret;
    ret.head<3>() = x.head<3>() + u + 0.01 * m.head<3>();
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
    ObserveT ret(observe_dim);
    ret(0) = ret(1) = ret(2) = 1;
    vector<Beam2D> cur_fov;
    truncate_beams(arm_helper->partitioned_fov, arm_helper->angle_to_segments(x.head<3>()), &cur_fov);
    if (approx_factor < 0) {
      if (inside(object_pos, cur_fov)) {
        ret(3) = ret(4) = 1;
      } else {
        ret(3) = ret(4) = 0;
      }
    } else {
      double dist = Geometry2D::sgndist(x.tail<2>(), cur_fov);
      double tol = 0.1;
      ret(3) = ret(4) = 1 - sigmoid(approx_factor * (dist + tol));
    }
    return ret;
  }

  ObserveT ArmObserveFunc::observation_masks(const StateT& x, double approx_factor) const {
    return observe_masks_from_object_position(x, x.tail<2>(), approx_factor);
  }

  ObserveT ArmObserveFunc::real_observation(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT ret;
    ret.head<3>() = x.head<3>() + 0.01 * n.head<3>();
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

  void ArmPlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
    assert(distmap != NULL);
    vector<Beam2D> cur_fov;
    vector<Segment> segs;
    if (state != NULL) {
      vector<Segment> segs = arm_helper->angle_to_segments(state->head<3>());
      truncate_beams(arm_helper->partitioned_fov, segs, &cur_fov);
    } else {
      cur_fov = arm_helper->partitioned_fov;
    }

    for (int j = 0; j < distmap->height(); ++j) {
      QRgb *line = (QRgb*) distmap->scanLine(j);
      for (int i = 0; i < distmap->width(); ++i) {
        double x = unscale_x(i),
               y = unscale_y(j);
        double dist = sgndist(Point(x, y), cur_fov);
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

  ArmOptPlotter::ArmOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   ArmPlotter(x_min, x_max, y_min, y_max, helper, parent),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void ArmOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    //if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
    //  // replot distmap
    //  distmap = QImage(width(), height(), QImage::Format_RGB32);
    //  if (states.size() > 0) {
    //    compute_distmap(&distmap, &states[0], arm_helper->belief_func->approx_factor);
    //  } else {
    //    compute_distmap(&distmap, NULL, arm_helper->belief_func->approx_factor);
    //  }
    //}
    if (states.size() > 0) {
      //painter.drawImage(0, 0, distmap);
      QPen robot_start_pen(QColor(255, 0, 0, 255), 4, Qt::SolidLine);
      QPen robot_middle_pen(QColor(255, 0, 0, 100), 4, Qt::SolidLine);
      QPen robot_end_pen(QColor(255, 0, 0, 255), 4, Qt::SolidLine);

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
    int T = 10;
    bool plotting = false;
    int method = 2;
    double noise_level = 1;
    double base_vec_array[] = {0, -1, 0};
    double start_vec_array[] = {PI/4+PI/16, PI/2, PI/4+PI/16, 0, 5};
    vector<double> base_vec(base_vec_array, end(base_vec_array));
    vector<double> start_vec(start_vec_array, end(start_vec_array));
    {
      Config config;
      config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter<int>("method", &method, "method"));
      config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }
    Vector5d start = toVectorXd(start_vec);
    Vector3d base_config = toVectorXd(base_vec);
    Matrix5d start_sigma = Matrix5d::Identity() * 0.01;
    start_sigma.bottomRightCorner<2, 2>() = Matrix2d::Identity() * 10;

    deque<Vector3d> initial_controls;
    for (int i = 0; i < T; ++i) {
      initial_controls.push_back(Vector3d::Zero());
    }

    Vector3d link_lengths; link_lengths << 4, 3, 2;

    Beam2D camera;
    camera.base = Point(0, 0);
    camera.a = Point(5, 10);
    camera.b = Point(-5, 10);

    ArmBSPPlannerPtr planner(new ArmBSPPlanner());

    planner->start = start;
    planner->base_config = base_config;
    planner->start_sigma = start_sigma;
    planner->link_lengths = link_lengths;
    planner->T = T;
    planner->noise_level = noise_level;
    planner->method = method;
    planner->controls = initial_controls;
    planner->camera = camera;
    planner->n_fov_parts = 10;
    planner->real_object_pos = Vector2d(0, 5);//2, 7);
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
      opt_callback = boost::bind(&ArmOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    //if (plotting) {
    //  emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
    //}

    planner->solve(opt_callback);
    StateT final_pos = getVec(planner->result, planner->helper->state_vars.row(planner->T));
    deque<Vector3d> _initial_controls;
    for (int i = 0; i < T; ++i) {
      _initial_controls.push_back((final_pos.head<3>() - start.head<3>()) / T);
    }
    planner->controls = planner->helper->initial_controls = _initial_controls;
    planner->simulate_executions(planner->T);
    planner->helper->real_object_pos = Vector2d(2, 7);//2, 7);
    if (plotting) {
      emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);//, planner->finished());
    }

    while (!planner->finished()) {
      planner->solve(opt_callback);
      planner->simulate_execution();
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);//, planner->finished());
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
