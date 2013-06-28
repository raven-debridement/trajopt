#include "bsp/bsp.hpp"
#include "car.hpp"
#include <QApplication>
#include <QtCore>

using namespace BSP;

namespace CarBSP {

  CarBSPProblemHelper::CarBSPProblemHelper() : BSPProblemHelper<CarBeliefFunc>() {
    input_dt = 0.25;
    carlen = 0.5;
    set_state_dim(4);
    set_sigma_dof(10);
    set_observe_dim(2);
    set_control_dim(2);

    double state_lbs_array[] = {-5, -3, -PI, 0};
    double state_ubs_array[] = {5, 3, PI, 2};
    double control_lbs_array[] = {-PI*0.25, -4};
    double control_ubs_array[] = {PI*0.25, 4};
    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));
    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 10);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
  }

  void CarBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      prob.addLinearConstraint(exprSub(AffExpr(state_vars.at(T, i)), goal(i)), EQ);
    }
  }

  void CarBSPProblemHelper::init_control_values(vector<ControlT>* output_init_controls) const {
    assert (output_init_controls != NULL);
    ControlT control_step = ControlT::Zero(control_dim);// = Control;
    //control_step.resize(control_dim);
    //control_step(0) = (goal(0) - start(0)) / T;
    //control_step(1) = (goal(1) - start(1)) / T;
    output_init_controls->clear();
    for (int i = 0; i < T; ++i) {
      output_init_controls->push_back(control_step);
    }
  }

  CarStateFunc::CarStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  CarStateFunc::CarStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  StateT CarStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT new_x(state_dim);
    /* Euler integration */
    //new_x(0) = x(0) + car_helper->input_dt * x(3) * cos(x(2));
    //new_x(1) = x(1) + car_helper->input_dt * x(3) * sin(x(2));
    //new_x(2) = x(2) + car_helper->input_dt * x(3) * tan(u(0)) / car_helper->carlen;
    //new_x(3) = x(3) + car_helper->input_dt * u(1);
    /* RK4 integration */
    StateT xtmp(state_dim), x1(state_dim), x2(state_dim), x3(state_dim), x4(state_dim);
    xtmp = x;
    x1(0) = xtmp(3) * cos(xtmp(2));
    x1(1) = xtmp(3) * sin(xtmp(2));
    x1(2) = xtmp(3) * tan(u(0)) / car_helper->carlen;
    x1(3) = u(1);
    xtmp = x + 0.5*car_helper->input_dt*x1;
    x2(0) = xtmp(3) * cos(xtmp(2));
    x2(1) = xtmp(3) * sin(xtmp(2));
    x2(2) = xtmp(3) * tan(u(0)) / car_helper->carlen;
    x2(3) = u(1);
    xtmp = x + 0.5*car_helper->input_dt*x2;
    x3(0) = xtmp(3) * cos(xtmp(2));
    x3(1) = xtmp(3) * sin(xtmp(2));
    x3(2) = xtmp(3) * tan(u(0)) / car_helper->carlen;
    x3(3) = u(1);
    xtmp = x + car_helper->input_dt*x3;
    x4(0) = xtmp(3) * cos(xtmp(2));
    x4(1) = xtmp(3) * sin(xtmp(2));
    x4(2) = xtmp(3) * tan(u(0)) / car_helper->carlen;
    x4(3) = u(1);
    
    new_x = x + car_helper->input_dt/6.0*(x1+2.0*(x2+x3)+x4);

    double umag = u.squaredNorm();
    return new_x + 0.01*umag*m;
  }

  CarObserveFunc::CarObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  CarObserveFunc::CarObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  ObserveT CarObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT out(observe_dim);
    out(0) = x(0) + 0.1 * n(0);
    out(1) = x(1) + 0.1 * n(1);
    return out;
  }

  CarBeliefFunc::CarBeliefFunc() : BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT>(), tol(0.1), alpha(0.5) {}

  CarBeliefFunc::CarBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
    tol(0.1), alpha(0.5), BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT>(helper, f, h), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  bool CarBeliefFunc::sgndist(const Vector2d& x, Vector2d* dists) const {
    Vector2d p1; p1 << 0, 2;
    Vector2d p2; p2 << 0, 0;
    (*dists)(0) = (x.head<2>() - p1).norm() - 0.5;
    (*dists)(1) = (x.head<2>() - p2).norm() - 0.5;
    return (*dists)(0) < 0 || (*dists)(1) < 0;
  }

  ObserveMatT CarBeliefFunc::compute_gamma(const StateT& x) const {
    Vector2d dists;
    sgndist(x.head<2>(), &dists);
    double gamma1 = 1. - (1./(1.+exp(-alpha*(dists(0)+tol))));
    double gamma2 = 1. - (1./(1.+exp(-alpha*(dists(1)+tol))));
    ObserveMatT gamma(observe_dim, observe_dim);
    gamma << gamma1, 0,
             0, gamma2;
    return gamma;
  }

  ObserveStateGradT CarBeliefFunc::sensor_constrained_observe_state_gradient(const ObserveStateGradT& H, const StateT& x) const {
    return compute_gamma(x) * H; 
  }

  VarianceT CarBeliefFunc::sensor_constrained_variance_reduction(const VarianceT& reduction, const StateT& x) const {
    ObserveMatT gamma = compute_gamma(x);
    VarianceT full_gamma = VarianceT::Identity(state_dim, state_dim);
    full_gamma.block(0, 0, observe_dim, observe_dim) = gamma;
    return full_gamma * reduction;
  }

  CarPlotter::CarPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)),
   old_alpha(-1), cur_alpha(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void CarPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (cur_alpha != old_alpha || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      for (int j = 0; j < height(); ++j) {
        QRgb *line = (QRgb*) distmap.scanLine(j);
        for (int i = 0; i < width(); ++i) {
          double x = unscale_x(i),
                 y = unscale_y(j);
          Vector2d dists;
          Vector2d pos; pos << x, y;
          car_helper->belief_func->sgndist(pos, &dists);
          double grayscale = fmax(1./(1. + exp(car_helper->belief_func->alpha*dists(0))),
                                  1./(1. + exp(car_helper->belief_func->alpha*dists(1))));
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
      draw_ellipse(states[i].head<2>(), sigmas[i].block<2,2>(0,0), painter, 0.5);
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

  void CarPlotter::update_plot_data(OptProb*, DblVec& xvec) {
    vector<VectorXd> new_states;
    vector<MatrixXd> new_sigmas;
    old_alpha = cur_alpha;
    cur_alpha = car_helper->belief_func->alpha;
    BeliefT cur_belief;
    car_helper->belief_func->compose_belief(car_helper->start, car_helper->start_sigma, &cur_belief);
    for (int i = 0; i <= T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      car_helper->belief_func->extract_state(cur_belief, &cur_state);
      car_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < T) cur_belief = car_helper->belief_func->call(cur_belief, (ControlT) getVec(xvec, car_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;
    this->repaint();
  }

  CarOptimizerTask::CarOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  CarOptimizerTask::CarOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void CarOptimizerTask::emit_plot_message(const boost::any& data) {
    BSPOptimizerTask::emit_plot_message(data);
  }

  void CarOptimizerTask::run() {
    int T = 20;

    bool plotting = false;

    double start_vec_array[] = {-5, 2, -PI*0.5, 0.2};
    double goal_vec_array[] = {-5, 0.0, 0.0, 0.0};

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

    Vector4d start = toVectorXd(start_vec);
    Vector4d goal = toVectorXd(goal_vec);
    Matrix4d start_sigma = Matrix4d::Identity(); //start_sigma << 1, 0, 0, 0,
                          //               0, 1, 0, 0,
                          //               0, 0, 1, 0,
                          //               0, 0, 0, 1;

    OptProbPtr prob(new OptProb());

    CarBSPProblemHelperPtr helper(new CarBSPProblemHelper());
    helper->start = start;
    helper->goal = goal;
    helper->start_sigma = start_sigma;
    helper->T = T;
    helper->configure_problem(*prob);

    BSPTrustRegionSQP opt(prob);
    opt.max_iter_ = 500;
    opt.merit_error_coeff_ = 500;
    opt.trust_shrink_ratio_ = 0.5;
    opt.trust_expand_ratio_ = 1.25;
    opt.min_trust_box_size_ = 1e-3;
    opt.min_approx_improve_ = 1e-2;
    opt.min_approx_improve_frac_ = 1e-4;
    opt.improve_ratio_threshold_ = 0.2;
    opt.trust_box_size_ = 1;

    helper->configure_optimizer(*prob, opt);

    boost::shared_ptr<CarPlotter> plotter;
    if (plotting) {
      double x_min = -7, x_max = 7, y_min = -5, y_max = 5;
      plotter.reset(create_plotter<CarPlotter>(x_min, x_max, y_min, y_max, helper));
      plotter->show();
      opt.addCallback(boost::bind(&CarOptimizerTask::emit_plot_message, boost::ref(this), _1, _2));
    }
    opt.optimize();
    if (!plotting) {
      emit finished_signal();
    }
  }
}

using namespace CarBSP;

int main(int argc, char *argv[]) {
  QApplication app(argc, argv);
  CarOptimizerTask* task = new CarOptimizerTask(argc, argv, &app);
  QTimer::singleShot(0, task, SLOT(run_slot()));
  QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
  return app.exec();
}
