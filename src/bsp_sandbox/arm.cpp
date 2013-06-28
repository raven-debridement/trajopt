#include "arm.hpp"
#include <QApplication>
#include <QtCore>

using namespace BSP;

namespace ArmBSP {

  class ArmBSPPlanner : public BSPPlanner<ArmBSPProblemHelper> {
  public:
    Vector2d goal_pos;

    virtual void initialize() {
      BSPPlanner<ArmBSPProblemHelper>::initialize(); 
      helper->goal_pos = goal_pos;
    }
  };
  
  typedef boost::shared_ptr<ArmBSPPlanner> ArmBSPPlannerPtr;

  ArmBSPProblemHelper::ArmBSPProblemHelper() : BSPProblemHelper<ArmBeliefFunc>() { initialize(); }

  void ArmBSPProblemHelper::initialize() {
    input_dt = 1.0;
    set_state_dim(3);
    set_sigma_dof(6);
    set_observe_dim(3);
    set_control_dim(3);
    set_state_bounds(DblVec(3, -PI), DblVec(3, PI));
    set_control_bounds(DblVec(3, -PI/4), DblVec(3, PI/4));
    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 40);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
    belief_constraints.clear();
  }

  TransT ArmBSPProblemHelper::angle_to_transform(const StateT& angle) const {
    TransT mat(2, 3);
    double t0 = base_config(2) + angle(0);
    double t1 = t0 + angle(1);
    double t2 = t1 + angle(2);
    mat << cos(t0), cos(t1), cos(t2),
           sin(t0), sin(t1), sin(t2);
    return mat;
  }

  Vector8d ArmBSPProblemHelper::angle_to_endpoints(const StateT& angle) const {
    Vector3d l1; l1 << link_length(0), 0, 0;
    Vector3d l2; l2 << link_length(0), link_length(1), 0;
    Vector3d l3; l3 << link_length(0), link_length(1), link_length(2);
    TransT mat = angle_to_transform(angle);
    Vector6d res;
    res.segment<2>(0) = base_config.head<2>();
    res.segment<2>(2) = base_config.head<2>() + mat * l1;
    res.segment<2>(4) = base_config.head<2>() + mat * l2;
    res.segment<2>(6) = base_config.head<2>() + mat * l3;
    return res;
  }

  Vector2d ArmBSPProblemHelper::angle_to_pos(const StateT& angle) const {
    Vector3d l3; l3 << link_length(0), link_length(1), link_length(2);
    TransT mat = angle_to_transform(angle);
    return mat * l3 + base_config.head<2>();
  }

  //StateT ArmBSPProblemHelper::pos_to_angle(const Vector2d& x) const {

  //}

  ArmStateFunc::ArmStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ArmStateFunc::ArmStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  StateT ArmStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    return x + u + 0.01*m;
  }

  ArmObserveFunc::ArmObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ArmObserveFunc::ArmObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {}

  ObserveT ArmObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return x + 0.1 * compute_inverse_gamma(x, -1) * n;
  }

  ObserveT ArmObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const {
    return x + 0.1 * compute_inverse_gamma(x, approx_factor) * n;
  }

  bool ArmObserveFunc::sgndist(const Vector2d& x, Vector2d* dists) const {
    Vector2d p1; p1 << 0, 2;
    Vector2d p2; p2 << 0, 0;
    (*dists)(0) = (x - p1).norm() - 0.5;
    (*dists)(1) = (x - p2).norm() - 0.5;
    return (*dists)(0) < 0 || (*dists)(1) < 0;
  }

  ObserveMatT ArmObserveFunc::compute_gamma(const StateT& x, double approx_factor) const {
    double tol = 0.1;
    Vector2d dists;
    sgndist(arm_helper->angle_to_pos(x), &dists);
    double gamma1, gamma2, gamma3;
    if (approx_factor < 0) {
      gamma1 = dists(0) <= 0 ? 1 : 0;
      gamma2 = dists(0) <= 0 ? 1 : 0;
      gamma3 = dists(1) <= 0 ? 1 : 0;
    } else {
      gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      gamma3 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
    }
    ObserveMatT gamma(observe_dim, observe_dim);
    gamma << gamma1, 0, 0,
             0, gamma2, 0,
             0, 0, gamma3;
    return gamma;
  }

  ObserveMatT ArmObserveFunc::compute_inverse_gamma(const StateT& x, double approx_factor) const {
    double tol = 0.1;
    Vector2d dists;
    sgndist(arm_helper->angle_to_pos(x), &dists);
    double minval = 1e-4;
    double invgamma1, invgamma2, invgamma3;

    if (approx_factor < 0) {
      invgamma1 = dists(0) <= 0 ? 1 : 1/minval;
      invgamma2 = dists(0) <= 0 ? 1 : 1/minval;
      invgamma3 = dists(1) <= 0 ? 1 : 1/minval;
    } else {
      double gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      double gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      double gamma3 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
      invgamma1 = 1. / (fmax(gamma1, minval));
      invgamma2 = 1. / (fmax(gamma2, minval));
      invgamma3 = 1. / (fmax(gamma3, minval));
    }

    ObserveMatT invgamma(observe_dim, observe_dim);
    invgamma << invgamma1, 0, 0,
                0, invgamma2, 0,
                0, 0, invgamma3;
    return invgamma;
  }

  ArmBeliefFunc::ArmBeliefFunc() : BeliefFunc<ArmStateFunc, ArmObserveFunc, BeliefT>() {
    set_approx_factor(0.5);
  }

  ArmBeliefFunc::ArmBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
    BeliefFunc<ArmStateFunc, ArmObserveFunc, BeliefT>(helper, f, h), arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)) {
    set_approx_factor(0.5);
  }

  ArmOptPlotter::ArmOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void ArmOptPlotter::keyPressEvent(QKeyEvent*) {

  }

  void ArmOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      for (int j = 0; j < height(); ++j) {
        QRgb *line = (QRgb*) distmap.scanLine(j);
        for (int i = 0; i < width(); ++i) {
          double x = unscale_x(i),
                 y = unscale_y(j);
          Vector2d dists;
          Vector2d pos; pos << x, y;
          arm_helper->belief_func->h->sgndist(pos, &dists);
          double grayscale = fmax(1./(1. + exp(arm_helper->belief_func->approx_factor*dists(0))),
                                  1./(1. + exp(arm_helper->belief_func->approx_factor*dists(1))));
          line[i] = qRgb(grayscale*255, grayscale*255, grayscale*255);
        }
      }
    }
    painter.drawImage(0, 0, distmap);
    QPen cvx_cov_pen(Qt::red, 2, Qt::SolidLine);
    QPen path_pen(Qt::red, 2, Qt::SolidLine);
    QPen pos_pen(Qt::red, 8, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(cvx_cov_pen);
    for (int i = 0; i < states.size(); ++i) {
      draw_ellipse(states[i], sigmas[i], painter, 0.5);
    }
    painter.setPen(path_pen);
    for (int i = 0; i < states.size() - 1; ++i) {
      draw_line(states[i](0), states[i](1), states[i+1](0), states[i+1](1), painter);
    }
    painter.setPen(pos_pen);
    for (int i = 0; i < states.size(); ++i) {
      draw_point(states[i](0), states[i](1), painter);
    }
  }

  void ArmOptPlotter::update_plot_data(void* data) {
    DblVec* xvec = (DblVec* ) data;
    vector<VectorXd> new_states;
    vector<MatrixXd> new_sigmas;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = arm_helper->belief_func->approx_factor;
    BeliefT cur_belief;
    arm_helper->belief_func->compose_belief(arm_helper->start, arm_helper->start_sigma, &cur_belief);
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
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   arm_helper(boost::static_pointer_cast<ArmBSPProblemHelper>(helper)),
   distmap(400, 400, QImage::Format_RGB32) {}

  void ArmSimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      for (int j = 0; j < height(); ++j) {
        QRgb *line = (QRgb*) distmap.scanLine(j);
        for (int i = 0; i < width(); ++i) {
          double x = unscale_x(i),
                 y = unscale_y(j);
          Vector2d dists;
          Vector2d pos; pos << x, y;
          arm_helper->belief_func->h->sgndist(pos, &dists);
          double grayscale = fmax(dists(0) <= 0 ? 1. : 0.,
                                  dists(1) <= 0 ? 1. : 0.);
          line[i] = qRgb(grayscale*255, grayscale*255, grayscale*255);
        }
      }
    }
    painter.drawImage(0, 0, distmap);
    if (simulated_positions.size() <= 0) {
      return;
    }
    QPen path_pen(Qt::blue, 2, Qt::SolidLine);
    QPen pos_pen(Qt::blue, 8, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(path_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i+1](0), simulated_positions[i+1](1), painter);
    }
    painter.setPen(pos_pen);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](0), simulated_positions[i](1), painter);
    }

    QPen opt_cvx_cov_pen(Qt::red, 2, Qt::SolidLine);
    QPen opt_path_pen(Qt::red, 2, Qt::SolidLine);
    QPen opt_pos_pen(Qt::red, 8, Qt::SolidLine);
    painter.setPen(opt_cvx_cov_pen);
    for (int i = 0; i < states.size(); ++i) {
      draw_ellipse(states[i], sigmas[i], painter, 0.5);
    }
    painter.setPen(opt_path_pen);
    for (int i = 0; i < states.size() - 1; ++i) {
      draw_line(states[i](0), states[i](1), states[i+1](0), states[i+1](1), painter);
    }
    painter.setPen(opt_pos_pen);
    for (int i = 0; i < states.size(); ++i) {
      draw_point(states[i](0), states[i](1), painter);
    }
  }

  void ArmSimulationPlotter::update_plot_data(void* data_x, void* data_sim) {//DblVec& xvec, vector<StateT>& new_simulated_positions) {
    cout << "start updating data" << endl;
    //vector<StateT>* new_simulated_positions = ;
    simulated_positions = *((vector<StateT>*) data_sim);
    //cout << "new simulated position size: " << new_simulated_positions.size() << endl;
    //for (int i = 0; i < new_simulated_positions->size(); ++i) {
    //  simulated_positions.push_back(new_simulated_positions->at(i));
    //}
    DblVec* xvec = (DblVec* ) data_x;
    vector<VectorXd> new_states;
    vector<MatrixXd> new_sigmas;
    BeliefT cur_belief;
    arm_helper->belief_func->compose_belief(arm_helper->start, arm_helper->start_sigma, &cur_belief);
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
    cout << "repainting" << endl;

    this->repaint();
    cout << "repainted" << endl;
  }

  ArmOptimizerTask::ArmOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  ArmOptimizerTask::ArmOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void ArmOptimizerTask::stage_plot_callback(boost::shared_ptr<ArmOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
  }

  void ArmOptimizerTask::run() {
    int T = 20;
    bool plotting = true;
    double base_vec_array[] = {-5, 2, 0};
    double start_vec_array[] = {PI/4, 0, 0};
    double goal_vec_array[] = {-5, 0};
    vector<double> base_vec(base_vec_array, end(base_vec_array));
    vector<double> start_vec(start_vec_array, end(start_vec_array));
    vector<double> goal_vec(goal_vec_array, end(goal_vec_array));
    {
      Config config;
      config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }
    Vector3d start = toVectorXd(start_vec);
    Vector2d goal_pos = toVectorXd(goal_vec);
    Matrix3d start_sigma = Matrix3d::Identity() * 0.1;

    deque<Vector3d> initial_controls;
    for (int i = 0; i < T; ++i) {
      initial_controls.push_back(Vector3d::Zero());
    }
    //deque<Vector2d> initial_controls;
    //Vector2d control_step = (goal - start) / T;
    //for (int i = 0; i < T; ++i) {
    //  initial_controls.push_back(control_step);
    //}

    //typedef boost::shared_ptr< BSPPlanner<ArmBSPProblemHelper> > BSPPlannerPtr;

    ArmBSPPlannerPtr planner(new ArmBSPPlanner());//<ArmBSPProblemHelper>());

    planner->start = start;
    planner->goal_pos = goal_pos;
    planner->start_sigma = start_sigma;
    planner->T = T;
    planner->controls = initial_controls;
    planner->initialize();

    boost::shared_ptr<ArmSimulationPlotter> sim_plotter;
    boost::shared_ptr<ArmOptPlotter> opt_plotter;
    if (plotting) {
      double x_min = -10, x_max = 10, y_min = -10, y_max = 10;
      sim_plotter.reset(create_plotter<ArmSimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      opt_plotter.reset(create_plotter<ArmOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      opt_plotter->show();
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    if (plotting) {
      opt_callback = boost::bind(&ArmOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    while (!planner->finished()) {
      if (plotting) {
        planner->solve(opt_callback);
      }
      planner->simulate_execution();
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
        //emit replot_signal(planner->result, planner->simulated_positions);
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