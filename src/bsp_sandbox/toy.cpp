#include "toy.hpp"
#include <QApplication>
#include <QtCore>

using namespace BSP;

namespace ToyBSP {

  ToyBSPProblemHelper::ToyBSPProblemHelper() : BSPProblemHelper<ToyBeliefFunc>() {
    input_dt = 1.0;
    set_state_dim(6);
    set_sigma_dof(21);
    set_observe_dim(6);
    set_control_dim(6);
    set_state_bounds(DblVec(6, -10), DblVec(6, 10));
    set_control_bounds(concat(DblVec(2, -0.8), DblVec(4, -0.4)), concat(DblVec(2, 0.8), DblVec(4, 0.4)));
    set_variance_cost(VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(VarianceT::Identity(state_dim, state_dim) * 10);
    set_control_cost(ControlCostT::Identity(control_dim, control_dim));
    belief_constraints.clear();
  }

  void ToyBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      // box constraint of buffer goaleps around goal position
      Var optvar = state_vars.at(T, i);
      prob.addLinearConstraint(exprSub(AffExpr(optvar), goal(i) ), EQ);
      //prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i)+goaleps) ), INEQ);
      //prob.addLinearConstraint(exprSub(exprSub(AffExpr(0), AffExpr(optvar)), (-goal(i)+goaleps) ), EQ);
    }
  }

  ToyStateFunc::ToyStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ToyStateFunc::ToyStateFunc(BSPProblemHelperBasePtr helper) :
    StateFunc<StateT, ControlT, StateNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  StateT ToyStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT new_x(state_dim);
    new_x(0) = x(0) + u(0)*toy_helper->input_dt + sqrt(0.01*u(0)*u(0)*toy_helper->input_dt + 0.001) * m(0);
    new_x(1) = x(1) + u(1)*toy_helper->input_dt + sqrt(0.01*u(1)*u(1)*toy_helper->input_dt + 0.001) * m(1);
    new_x.tail<4>() = x.tail<4>() + u.tail<4>() + 0.01 * m.tail<4>();
    return new_x;
  }

  ToyObserveFunc::ToyObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ToyObserveFunc::ToyObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return x + 0.1 * compute_inverse_gamma(x, -1) * n;
  }

  ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const {
    return x + 0.1 * compute_inverse_gamma(x, approx_factor) * n;
  }

  bool ToyObserveFunc::sgndist(const StateT& x, Vector2d* dists) const {
    (*dists)(0) = (x.head<2>() - x.middleRows<2>(2)).norm() - 0.5;
    (*dists)(1) = (x.head<2>() - x.tail<2>()).norm() - 0.5;
    return (*dists)(0) < 0 || (*dists)(1) < 0;
  }

  ObserveMatT ToyObserveFunc::compute_gamma(const StateT& x, double approx_factor) const {
    double tol = 0.2;
    Vector2d dists;
    sgndist(x, &dists);
    double gamma1, gamma2;
    if (approx_factor < 0) {
      gamma1 = dists(0) <= 0 ? 1 : 0;
      gamma2 = dists(1) <= 0 ? 1 : 0;
    } else {
      gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
    }
    ObserveMatT gamma(observe_dim, observe_dim);
    gamma << gamma1, 0,      0, 0, 0, 0,
             0,      gamma2, 0, 0, 0, 0,
             0,      0,      1, 0, 0, 0,
             0,      0,      0, 1, 0, 0,
             0,      0,      0, 0, 1, 0,
             0,      0,      0, 0, 0, 1;
    return gamma;
  }

  ObserveMatT ToyObserveFunc::compute_inverse_gamma(const StateT& x, double approx_factor) const {
    double tol = 0.2;
    Vector2d dists;
    sgndist(x, &dists);
    double minval = 1e-4;
    double invgamma1, invgamma2;

    if (approx_factor < 0) {
      invgamma1 = dists(0) <= 0 ? 1 : 1/minval;
      invgamma2 = dists(1) <= 0 ? 1 : 1/minval;
    } else {
      double gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      double gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
      invgamma1 = 1. / (fmax(gamma1, minval));
      invgamma2 = 1. / (fmax(gamma2, minval));
    }

    ObserveMatT invgamma(observe_dim, observe_dim);
    invgamma << invgamma1, 0,         0, 0, 0, 0,
                0,         invgamma2, 0, 0, 0, 0,
                0,         0,         1, 0, 0, 0,
                0,         0,         0, 1, 0, 0,
                0,         0,         0, 0, 1, 0,
                0,         0,         0, 0, 0, 1;
    return invgamma;
  }

  ToyBeliefFunc::ToyBeliefFunc() : BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>() {
    set_approx_factor(0.5);
  }

  ToyBeliefFunc::ToyBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
    BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>(helper, f, h), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {
    set_approx_factor(0.5);
  }

  ToyPlotter::ToyPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  void ToyPlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
    assert(distmap != NULL);
    if (state == NULL) {
      for (int j = 0; j < distmap->height(); ++j) {
        QRgb *line = (QRgb*) distmap->scanLine(j);
        for (int i = 0; i < distmap->width(); ++i) {
          line[i] = qRgb(0, 0, 0);
        }
      }
      return;
    }
    for (int j = 0; j < distmap->height(); ++j) {
      QRgb *line = (QRgb*) distmap->scanLine(j);
      for (int i = 0; i < distmap->width(); ++i) {
        double x = unscale_x(i),
               y = unscale_y(j);
        StateT pos;
        pos.head<2>() = Vector2d(x, y);
        pos.tail<4>() = state->tail<4>();
        Vector2d dists;
        toy_helper->belief_func->h->sgndist(pos, &dists);
        double grayscale;
        if (approx_factor > 0) {
          grayscale = fmax(1./(1. + exp(approx_factor*dists(0))),
                           1./(1. + exp(approx_factor*dists(1))));
        } else {
          grayscale = dists(0) <= 0 || dists(1) <= 0 ? 1 : 0;
        }
        line[i] = qRgb(grayscale*255, grayscale*255, grayscale*255);
      }
    }
  }

  ToyOptPlotter::ToyOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   ToyPlotter(x_min, x_max, y_min, y_max, helper, parent),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void ToyOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      if (states_opt.size() > 0) {
        compute_distmap(&distmap, &states_opt[0], toy_helper->belief_func->approx_factor);
      } else {
        compute_distmap(&distmap, NULL, toy_helper->belief_func->approx_factor);
      }
    }
    //painter.drawImage(0, 0, distmap);
    QPen cvx_cov_pen(Qt::red, 3, Qt::SolidLine);
    QPen path_pen(Qt::red, 3, Qt::SolidLine);
    QPen pos_pen(Qt::red, 8, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(cvx_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].head<2>(), sigmas_opt[i].topLeftCorner<2, 2>(), painter, 0.5);
    }
    painter.setPen(path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](0), states_opt[i](1), states_opt[i+1](0), states_opt[i+1](1), painter);
    }
    painter.setPen(pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](0), states_opt[i](1), painter);
    }

    QPen sensor_cov_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].middleRows<2>(2), sigmas_opt[i].block<2, 2>(2, 2), painter, 0.5);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](2), states_opt[i](3), states_opt[i+1](2), states_opt[i+1](3), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](2), states_opt[i](3), painter);
    }

    painter.setPen(sensor_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].middleRows<2>(4), sigmas_opt[i].block<2, 2>(4, 4), painter, 0.5);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](4), states_opt[i](5), states_opt[i+1](4), states_opt[i+1](5), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](4), states_opt[i](5), painter);
    }
    // draw beliefs computed using belief dynamics
    //QPen cvx_cov_pen2(Qt::blue, 3, Qt::SolidLine);
    //painter.setPen(cvx_cov_pen2);
    //for (int i = 0; i < states_actual.size(); ++i) {
    //  draw_ellipse(states_actual[i], sigmas_actual[i], painter, 0.5);
    //}
    //QPen pos_pen2(Qt::blue, 8, Qt::SolidLine);
    //painter.setPen(pos_pen2);
    //for (int i = 0; i < states_actual.size(); ++i) {
    //  draw_point(states_actual[i](0), states_actual[i](1), painter);
    //}
  }

  void ToyOptPlotter::update_plot_data(void* data) {
    DblVec* xvec = (DblVec* ) data;
    vector<StateT> new_states_opt, new_states_actual;
    vector<VarianceT> new_sigmas_opt, new_sigmas_actual;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = toy_helper->belief_func->approx_factor;
    toy_helper->belief_func->approx_factor = 1000;
    BeliefT cur_belief_opt, cur_belief_actual;
    toy_helper->belief_func->compose_belief(toy_helper->start, matrix_sqrt(toy_helper->start_sigma), &cur_belief_actual);
    int T = toy_helper->get_T();

    for (int i = 0; i <= T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;

      cur_belief_opt = (BeliefT) getVec(*xvec, toy_helper->belief_vars.row(i));
      toy_helper->belief_func->extract_state(cur_belief_opt, &cur_state);
      toy_helper->belief_func->extract_sigma(cur_belief_opt, &cur_sigma);
      new_states_opt.push_back(cur_state);
      new_sigmas_opt.push_back(cur_sigma);

      toy_helper->belief_func->extract_state(cur_belief_actual, &cur_state);
      toy_helper->belief_func->extract_sigma(cur_belief_actual, &cur_sigma);
      new_states_actual.push_back(cur_state);
      new_sigmas_actual.push_back(cur_sigma);
      if (i < T) {
        cur_belief_actual = toy_helper->belief_func->call(cur_belief_actual, (ControlT) getVec(*xvec, toy_helper->control_vars.row(i)));
      }
    }

    states_opt = new_states_opt; states_actual = new_states_actual;
    sigmas_opt = new_sigmas_opt; sigmas_actual = new_sigmas_actual;
    toy_helper->belief_func->approx_factor = cur_approx_factor;
    this->repaint();
  }

  ToySimulationPlotter::ToySimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   ToyPlotter(x_min, x_max, y_min, y_max, helper, parent),
   distmap(400, 400, QImage::Format_RGB32) {}

  void ToySimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      if (simulated_positions.size() > 0) {
        compute_distmap(&distmap, &simulated_positions[0], -1);
      } else {
        compute_distmap(&distmap, NULL, -1);
      }
    }
    //painter.drawImage(0, 0, distmap);
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

    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](2), simulated_positions[i](3), simulated_positions[i+1](2), simulated_positions[i+1](3), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](2), simulated_positions[i](3), painter);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](4), simulated_positions[i](5), simulated_positions[i+1](4), simulated_positions[i+1](5), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](4), simulated_positions[i](5), painter);
    }
    //QPen opt_cvx_cov_pen(Qt::red, 3, Qt::SolidLine);
    //QPen opt_path_pen(Qt::red, 3, Qt::SolidLine);
    //QPen opt_pos_pen(Qt::red, 8, Qt::SolidLine);
    //painter.setPen(opt_cvx_cov_pen);
    //for (int i = 0; i < states.size(); ++i) {
    //  draw_ellipse(states[i], sigmas[i], painter, 0.5);
    //}
    //painter.setPen(opt_path_pen);
    //for (int i = 0; i < states.size() - 1; ++i) {
    //  draw_line(states[i](0), states[i](1), states[i+1](0), states[i+1](1), painter);
    //}
    //painter.setPen(opt_pos_pen);
    //for (int i = 0; i < states.size(); ++i) {
    //  draw_point(states[i](0), states[i](1), painter);
    //}
  }

  void ToySimulationPlotter::update_plot_data(void* data_x, void* data_sim) {
    simulated_positions = *((vector<StateT>*) data_sim);
    DblVec* xvec = (DblVec* ) data_x;
    vector<StateT> new_states;
    vector<VarianceT> new_sigmas;
    BeliefT cur_belief;
    toy_helper->belief_func->compose_belief(toy_helper->start, matrix_sqrt(toy_helper->start_sigma), &cur_belief);
    for (int i = 0; i <= toy_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      toy_helper->belief_func->extract_state(cur_belief, &cur_state);
      toy_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < toy_helper->T) cur_belief = toy_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, toy_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;

    this->repaint();
  }

  ToyOptimizerTask::ToyOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  ToyOptimizerTask::ToyOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void ToyOptimizerTask::stage_plot_callback(boost::shared_ptr<ToyOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
  }

  void ToyOptimizerTask::run() {
    int T = 20;
    bool plotting = true;
    double start_vec_array[] = {-5, 2, 0, 2, 0, 0};
    double goal_vec_array[] = {-5, 0, 0, 0, 0, 0};
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
    Vector6d start = toVectorXd(start_vec);
    Vector6d goal = toVectorXd(goal_vec);
    Matrix6d start_sigma = Matrix6d::Identity();

    deque<Vector6d> initial_controls;
    Vector6d control_step = (goal - start) / T;
    control_step.tail<4>() = Vector4d::Zero();
    for (int i = 0; i < T; ++i) {
      initial_controls.push_back(control_step);
    }

    typedef boost::shared_ptr< BSPPlanner<ToyBSPProblemHelper> > BSPPlannerPtr;

    BSPPlannerPtr planner(new BSPPlanner<ToyBSPProblemHelper>());

    planner->start = start;
    planner->goal = goal;
    planner->start_sigma = start_sigma;
    planner->T = T;
    planner->controls = initial_controls;
    planner->initialize();

    boost::shared_ptr<ToySimulationPlotter> sim_plotter;
    //boost::shared_ptr<ToyOptPlotter> opt_plotter;
    if (plotting) {
      double x_min = -10, x_max = 10, y_min = -10, y_max = 10;
      sim_plotter.reset(create_plotter<ToySimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      //opt_plotter.reset(create_plotter<ToyOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      //opt_plotter->show();
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    //if (plotting) {
    //  opt_callback = boost::bind(&ToyOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    //}

    while (!planner->finished()) {
      planner->solve(opt_callback);
      planner->simulate_execution();
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
      }
    }
    emit finished_signal();
  }
}

using namespace ToyBSP;

int main(int argc, char *argv[]) {
	QApplication app(argc, argv);
	ToyOptimizerTask* task = new ToyOptimizerTask(argc, argv, &app);
	QTimer::singleShot(0, task, SLOT(run_slot()));
	QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
	return app.exec();
}
