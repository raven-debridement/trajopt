#include "bsp/bsp.hpp"
#include "toy.hpp"
#include <deque>
#include <QApplication>
#include <QtCore>
#include <QPolygonF>


using namespace BSP;

namespace ToyBSP {

  ToyBSPPlanner::ToyBSPPlanner() : BSPPlanner<ToyBSPProblemHelper>() {
    this->n_alpha_iterations = 5;
  }

  void ToyBSPPlanner::custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state) {
    assert (state != NULL);
    assert (sigma != NULL);
    Vector2d new_state;
    Matrix2d new_sigma;

    if (actual_state(0) > 0) {
      truncate_gaussian(Vector2d(-1, 0), 0, *state, *sigma, &new_state, &new_sigma);
    } else {
      truncate_gaussian(Vector2d(1, 0), 0, *state, *sigma, &new_state, &new_sigma);
    }

    //std::cout << "state: " << (*state)(0) << " " << (*state)(1) << " sigma: " << (*sigma)(0,0) << " " << (*sigma)(0,1) << " " << (*sigma)(1,0) << " " << (*sigma)(1,1) << endl;

    //if ((*state)(0) > 0 && actual_state(0) < 0) {
    //	truncate_gaussian(Vector2d(-1, 0), 0, *state, *sigma, &new_state, &new_sigma);
    //} else if ((*state)(0) < 0 && actual_state(0) > 0){
    //	truncate_gaussian(Vector2d(1, 0), 0, *state, *sigma, &new_state, &new_sigma);
    //}
    //cout << "new_state: " << new_state(0) << " " << new_state(1) << endl;

    //*state = new_state;
    //*sigma = new_sigma;
  }

  void ToyBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt) {
    opt.max_iter_                   = 250;
		opt.merit_error_coeff_          = 100;
		opt.merit_coeff_increase_ratio_ = 10;
		opt.max_merit_coeff_increases_  = 3;
		opt.trust_shrink_ratio_         = 0.9;
		opt.trust_expand_ratio_         = 1.1;
		opt.min_trust_box_size_         = 1e-4;
		opt.min_approx_improve_         = 0.01;
		opt.min_approx_improve_frac_    = 1e-4;
		opt.improve_ratio_threshold_    = 0.1;
		opt.trust_box_size_             = 1;
		opt.cnt_tolerance_              = 1e-08;

  }

  ToyBSPProblemHelper::ToyBSPProblemHelper() : BSPProblemHelper<ToyBeliefFunc>() {
    input_dt = 1;
    goaleps = 0.25;

    set_state_dim(2);
    set_sigma_dof(3);
    set_observe_dim(2);
    set_control_dim(2);

    double state_lbs_array[] = {-7, -5};
    double state_ubs_array[] = {3, 5};

    double control_lbs_array[] = {-0.9, -0.9};
    double control_ubs_array[] = {0.9, 0.9};

    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));

    VarianceT variance_cost = VarianceT::Identity(state_dim, state_dim);
    set_variance_cost(variance_cost);
    set_final_variance_cost(variance_cost * 10);

    ControlCostT control_cost = ControlCostT::Identity(control_dim, control_dim);
    set_control_cost(control_cost);
  }

  void ToyBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      // box constraint of buffer goaleps around goal position
      Var optvar = state_vars.at(T, i);
      prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i))), EQ);
      //prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i)+goaleps) ), INEQ);
      //prob.addLinearConstraint(exprSub(exprSub(AffExpr(0), AffExpr(optvar)), (-goal(i)+goaleps) ), INEQ);
    }
  }

  ToyStateFunc::ToyStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  ToyStateFunc::ToyStateFunc(BSPProblemHelperBasePtr helper) :
                                      StateFunc<StateT, ControlT, StateNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  StateT ToyStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    //return x + u * toy_helper->input_dt + 0.05 * sqrt((u.transpose()*u).trace())* m;
    return x + u * toy_helper->input_dt + 0.05 * m;
  }

  ToyObserveFunc::ToyObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  ToyObserveFunc::ToyObserveFunc(BSPProblemHelperBasePtr helper) :
                        ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  ObserveT ToyObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return x + 0.01 * n;
  }

  ObserveT ToyObserveFunc::observation_masks(const StateT& x, double approx_factor) const {
    Vector1d dists;
    sgndist(x, &dists);
    double tol = 0.1;
    ObserveT ret(observe_dim);
    if (approx_factor < 0) {
      ret(0) = ret(1) = dists(0) <= 0 ? 1 : 0;
    } else {
      ret(0) = ret(1) = 1. - sigmoid(approx_factor * (dists(0) + tol));
    }
  }

  ToyBeliefFunc::ToyBeliefFunc() : BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>() {
    this->approx_factor = 1;
  }

  ToyBeliefFunc::ToyBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
                       BeliefFunc<ToyStateFunc, ToyObserveFunc, BeliefT>(helper, f, h), toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {
    this->approx_factor = 1;
  }

  ToyPlotter::ToyPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
                       BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
                       ProblemState(helper),
                       toy_helper(boost::static_pointer_cast<ToyBSPProblemHelper>(helper)) {}

  void ToyPlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
    assert(distmap != NULL);
    for (int j = 0; j < distmap->height(); ++j) {
      QRgb *line = (QRgb*) distmap->scanLine(j);
      for (int i = 0; i < distmap->width(); ++i) {
        if (state == NULL) {
          line[i] = qRgb(0, 0, 0);
        } else {
          double x = unscale_x(i),
              y = unscale_y(j);
          Vector2d pos; pos << x, y;
          Vector1d dists;
          toy_helper->belief_func->h->sgndist((StateT) pos, &dists);
          double grayscale;
          if (approx_factor > 0) {
            grayscale = 1./(1. + exp(approx_factor*dists(0)));
          } else {
            grayscale = dists(0) <= 0 ? 1 : 0;
          }
          line[i] = qRgb(grayscale*255, grayscale*255, grayscale*255);
        }
      }
    }
  }

  ToyOptPlotter::ToyOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
                       ToyPlotter(x_min, x_max, y_min, y_max, helper, parent),
                       old_approx_factor(-1), cur_approx_factor(-1), fileinp(false), distmap(400, 400, QImage::Format_RGB32) {}

  void ToyOptPlotter::paintEvent(QPaintEvent* ) {
    //cout << "cur_approx_factor: " << cur_approx_factor << endl;
    QPainter painter(this);
    if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      //cout << "approx_factor inside paintEvent: " << toy_helper->belief_func->approx_factor << endl;
      if ((int)states_opt.size() > 0) {
        compute_distmap(&distmap, &states_opt[0], toy_helper->belief_func->approx_factor);
      } else {
        compute_distmap(&distmap, NULL, toy_helper->belief_func->approx_factor);
      }
    }

    double ellipseScaling = 1.0;
    painter.drawImage(0, 0, distmap);
    QPen cvx_cov_pen(Qt::red, 3, Qt::SolidLine);
    //QPen path_pen(Qt::blue, 3, Qt::SolidLine);
    QPen path_pen(Qt::red, 3, Qt::SolidLine);
    QPen openlooppath_pen(Qt::green, 3, Qt::SolidLine);
    QPen pos_pen(Qt::red, 4, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);


    // draw goal
    QPen goal_pen(QColor(192,192,192,100), 2, Qt::SolidLine);
    painter.setPen(goal_pen);
    QBrush goal_brush(QColor(192, 192, 192, 100));
    QBrush prev_brush = painter.brush();
    painter.setBrush(goal_brush);
    Vector2d g = toy_helper->goal.head<2>();
    draw_ellipse(g, Matrix2d::Identity()*toy_helper->goaleps, painter, ellipseScaling);
    painter.setBrush(prev_brush);

    painter.setPen(cvx_cov_pen);
    //for (int i = 0; i < (int)states_opt.size(); ++i) {
    //	draw_ellipse(states_opt[i].head<2>(), sigmas_opt[i].topLeftCorner(2,2), painter, ellipseScaling);
    //}

    painter.setPen(path_pen);
    for (int i = 0; i < (int)states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](0), states_opt[i](1), states_opt[i+1](0), states_opt[i+1](1), painter);
    }

    /*
    painter.setPen(openlooppath_pen);
    for (int i = 0; i < 20; ++i) {
      for(int j = 0; j < 20; ++j) {
        draw_line(openlooppts[i*21+j](0), openlooppts[i*21+j](1), openlooppts[i*21+j+1](0), openlooppts[i*21+j+1](1),painter);
      }
    }
    */

    //painter.setPen(pos_pen);
    //for (int i = 0; i < (int)states_opt.size(); ++i) {
    // draw_point(states_opt[i](0), states_opt[i](1), painter);
    //}//truncate_gaussian(Vector2d(-1, 0), 0, *state, *sigma, &new_state, &new_sigma);


    // draw beliefs computed using belief dynamics
    QPen cvx_cov_pen2(QColor(255,215,0), 3, Qt::SolidLine);
    painter.setPen(cvx_cov_pen2);
    draw_ellipse(states_actual[0].head<2>(), sigmas_actual[0].topLeftCorner(2,2), painter, ellipseScaling);
    draw_ellipse(states_actual[(int)states_actual.size()-1].head<2>(), sigmas_actual[(int)states_actual.size()-1].topLeftCorner(2,2), painter, ellipseScaling);

    //for (int i = 0; i < (int)states_actual.size(); ++i) {
    //	draw_ellipse(states_actual[i].head<2>(), sigmas_actual[i].topLeftCorner(2,2), painter, ellipseScaling);
    //}
    //for (int i = 0; i < (int)states_01.size(); ++i) {
    //	draw_ellipse(states_01[i].head<2>(), sigmas_01[i].topLeftCorner(2,2), painter, ellipseScaling);
    //}

    //QPen pos_pen2(Qt::red, 8, Qt::SolidLine);
    QPen pos_pen2(QColor(255,215,0), 8, Qt::SolidLine);
    painter.setPen(pos_pen2);
    //for (int i = 0; i < (int)states_actual.size(); ++i) {
    //	draw_point(states_actual[i](0), states_actual[i](1), painter);
    //}
    for (int i = 0; i < (int)states_01.size(); ++i) {
      draw_point(states_01[i](0), states_01[i](1), painter);
    }

    QPen real_pen(QColor(0,255,0), 3, Qt::SolidLine);
    painter.setPen(real_pen);
    if ((int)simplotptr->simulated_positions.size() <= 0) {
      return;
    }
    for (int i = 0; i < (int)simplotptr->simulated_positions.size() - 1; ++i) {
      draw_line(simplotptr->simulated_positions[i](0), simplotptr->simulated_positions[i](1), simplotptr->simulated_positions[i+1](0), simplotptr->simulated_positions[i+1](1), painter);
    }
  }

  void ToyOptPlotter::update_plot_data(void* data)
  {
    vector<StateT> new_states_opt, new_states_actual, new_states_01;
    vector<VarianceT> new_sigmas_opt, new_sigmas_actual, new_sigmas_01;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = toy_helper->belief_func->approx_factor;
    //cout << "approx_factor inside update_plot_data: " << toy_helper->belief_func->approx_factor << endl;
    DblVec* xvec = (DblVec* ) data;
    BeliefT cur_belief_opt, cur_belief_actual, cur_belief_01;
    toy_helper->belief_func->compose_belief(toy_helper->start, matrix_sqrt(toy_helper->start_sigma), &cur_belief_actual);
    toy_helper->belief_func->compose_belief(toy_helper->start, matrix_sqrt(toy_helper->start_sigma), &cur_belief_01);
    int T = toy_helper->get_T();
    states_waypoints.clear();

    //for (int i = 0; i < T; ++i) {
    //	ControlT u = (ControlT) getVec(*xvec, toy_helper->control_vars.row(i));
    //	cout << u(0) << " " << u(1) << endl;
    //}
    //cout << endl;

    /*
    if (!fileinp) {
      ifstream ifptr("toy-openloop-paths.txt", ios::in);
      StateT xvec;
      for(int i = 0; i < 440; ++i) {
        ifptr >> xvec(0) >> xvec(1);
        openlooppts.push_back(xvec);
      }
      cout << xvec(0) << " " << xvec(1) << endl;
      ifptr.close();
      fileinp = true;
    }
    */

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

      toy_helper->belief_func->extract_state(cur_belief_01, &cur_state);
      toy_helper->belief_func->extract_sigma(cur_belief_01, &cur_sigma);
      new_states_01.push_back(cur_state);
      new_sigmas_01.push_back(cur_sigma);

      if (i < T) {
        cur_belief_actual = toy_helper->belief_func->call(cur_belief_actual, (ControlT) getVec(*xvec, toy_helper->control_vars.row(i)));

        toy_helper->belief_func->approx_factor = -1;
        cur_belief_01 = toy_helper->belief_func->call(cur_belief_01, (ControlT) getVec(*xvec, toy_helper->control_vars.row(i)));
        toy_helper->belief_func->approx_factor = cur_approx_factor;
      }
    }
    states_opt = new_states_opt; states_actual = new_states_actual; states_01 = new_states_01;
    sigmas_opt = new_sigmas_opt; sigmas_actual = new_sigmas_actual; sigmas_01 = new_sigmas_01;

    toy_helper->belief_func->approx_factor = cur_approx_factor;
    this->repaint();
  }

  ToySimulationPlotter::ToySimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
                       ToyPlotter(x_min, x_max, y_min, y_max, helper, parent),
                       distmap(400, 400, QImage::Format_RGB32) {}

  void ToySimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    distmap = QImage(width(), height(), QImage::Format_RGB32);
    if (simulated_positions.size() <= 0) {
      return;
    }
    compute_distmap(&distmap, &simulated_positions.back(), -1);
    painter.drawImage(0, 0, distmap);
    QPen sim_prev_pos_pen(QColor(255, 0, 0, 100), 4, Qt::SolidLine);
    QPen sim_cur_pos_pen(QColor(255, 0, 0, 255), 4, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(sim_cur_pos_pen);
    for (int i = 0; i < (int)simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i+1](0), simulated_positions[i+1](1), painter);
    }
    // draw goal
    QPen goal_pen(QColor(192,192,192,100), 2, Qt::SolidLine);
    painter.setPen(goal_pen);
    QBrush goal_brush(QColor(192, 192, 192, 100));
    QBrush prev_brush = painter.brush();
    painter.setBrush(goal_brush);
    Vector2d g = toy_helper->goal.head<2>();
    draw_ellipse(g, Matrix2d::Identity()*toy_helper->goaleps, painter, 1);
    painter.setBrush(prev_brush);
  }

  void ToySimulationPlotter::update_plot_data(void* data_x, void* data_sim) {
    simulated_positions = *((vector<StateT>*) data_sim);
    DblVec* xvec = (DblVec* ) data_x;
    vector<StateT> new_states;
    vector<VarianceT> new_sigmas;
    BeliefT cur_belief;
    toy_helper->belief_func->compose_belief(toy_helper->start, matrix_sqrt(toy_helper->start_sigma), &cur_belief);
    waypoints.clear();

    for (int i = 0; i <= toy_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      toy_helper->belief_func->extract_state(cur_belief, &cur_state);
      toy_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      cur_belief = toy_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, toy_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;

    this->repaint();
  }

  ToyOptimizerTask::ToyOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  ToyOptimizerTask::ToyOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void ToyOptimizerTask::stage_plot_callback(boost::shared_ptr<ToyOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
    //wait_to_proceed(boost::bind(&ToyOptPlotter::update_plot_data, plotter, &x));
  }


  void ToyOptimizerTask::run() {
    bool plotting = true;

    int method = 2;
    double noise_level =1;

    double start_vec_array[] = {-5, 2};
    double goal_vec_array[] = {-5, -2};

    vector<double> start_vec(start_vec_array, end(start_vec_array));
    vector<double> goal_vec(goal_vec_array, end(goal_vec_array));
    int seed;
    {
      Config config;
      config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
      config.add(new Parameter<int>("method", &method, "method"));
      config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
      config.add(new Parameter<int>("seed", &seed, "seed"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }
    srand(static_cast<unsigned int>(std::time(0)));
    //srand(seed);

    Vector2d start = toVectorXd(start_vec);
    Vector2d goal = toVectorXd(goal_vec);
    Matrix2d start_sigma = Matrix2d::Identity();

    int T = 20;
    deque<Vector2d> initial_controls;
    //ifstream ifptr("toy-opt-init.txt", ios::in);
    //int nu;
    //ifptr >> nu;
    ControlT uvec;
    for (int i = 0; i < T; ++i) {
      //initial_controls.push_back(Vector2d::Zero());
      initial_controls.push_back((goal-start)/(double)T);
      //ifptr >> uvec(0) >> uvec(1);
      //initial_controls.push_back(uvec);
    }
    //if (ifptr.is_open()) ifptr.close();

    ToyBSPPlannerPtr planner(new ToyBSPPlanner());
    planner->start = start;
    planner->goal = goal;
    planner->start_sigma = start_sigma;
    planner->method = method;
    planner->T = T;
    planner->controls = initial_controls;
    planner->noise_level = noise_level;
    planner->initialize();

    boost::shared_ptr<ToySimulationPlotter> sim_plotter;
    boost::shared_ptr<ToyOptPlotter> opt_plotter;

    if (plotting) {
      double x_min = -8.4, x_max = 3.9, y_min = -3.1, y_max = 3.1;
      sim_plotter.reset(create_plotter<ToySimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      //if (method != 0) {
      opt_plotter.reset(create_plotter<ToyOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      opt_plotter->show();
      //}
      opt_plotter->simplotptr = sim_plotter;
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    //if (plotting && method != 0) {
    if (plotting) {
      opt_callback = boost::bind(&ToyOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    // test truncated gaussians
    /*
    StateT state;
    state(0) = -4.36543; state(1) = 1.80434;
    VarianceT sigma;
    sigma(0,0) = 1.0024; sigma(1,0) = 0; sigma(0,1) = 0; sigma(1,1) = 1.0024;

    Vector2d new_state;
    Matrix2d new_sigma;
    truncate_gaussian(Vector2d(-1, 0), 0, state, sigma, &new_state, &new_sigma);

    cout << "state: " << state(0) << " " << state(1) << endl;
    cout << "sigma: " << sigma(0,0) << " " << sigma(0,1) << " " << sigma(1,0) << " " << sigma(1,1) << endl;

    cout << "new_state: " << new_state(0) << " " << new_state(1) << endl;
    cout << "new_sigma: " << new_sigma(0,0) << " " << new_sigma(0,1) << " " << new_sigma(1,0) << " " << new_sigma(1,1) << endl;

    cout << "PAUSED: testing truncated Gaussians" << endl;
    int num;
    cin >> num;
    */

    //cout << "start solving" << endl;

    /*
    int success = 0;
    for (int i = 0; i < 1000; ++i) {
      planner->helper->T = 20;
      initial_controls.clear();
      ifstream ifptr("toy-opt-init.txt", ios::in);
      int nu;
      ifptr >> nu;
      ControlT uvec;
      for (int i = 0; i < T; ++i) {
        ifptr >> uvec(0) >> uvec(1);
        initial_controls.push_back(uvec);
      }
      if (ifptr.is_open()) ifptr.close();
      planner->controls = initial_controls;

      planner->start = start;
      planner->goal = goal;
      planner->start_sigma = start_sigma;
      planner->method = method;
      planner->T = T;
      planner->noise_level = noise_level;
      planner->initialize();

      planner->simulate_executions(planner->helper->T);
      cout << "final position: " << planner->simulated_positions.back().head<2>().transpose() << endl;
      cout << "final estimated position: " << planner->start.head<2>().transpose() << endl;
      cout << "final covariance: " << endl << planner->start_sigma.topLeftCorner<2, 2>() << endl;
      cout << "goal position: " << planner->goal.head<2>().transpose() << endl;

      Vector2d fpos = planner->simulated_positions.back().head<2>();
      Vector2d goalpos = planner->goal.head<2>();
      double dist = (fpos - goalpos).norm();
      cout << "dist: " << dist << endl;
      cout << "goaleps: " << planner->helper->goaleps << endl;
      if (dist < planner->helper->goaleps) {
        cout << "success" << endl;
        ++success;
      }
      cout << "Successful runs: " << success << "/1000" << endl;
    }

    if (plotting) {
      emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
    }
    */

    int success =0;
    for (int i = 0; i < 20; ++i) {
      planner->helper->T = 20;
      planner->start = start;
      planner->goal = goal;
      ControlT uvec;
      initial_controls.clear();
      for (int i = 0; i < T; ++i) {
        initial_controls.push_back((goal-start)/(double)T);
      }
      planner->controls = initial_controls;
      planner->start_sigma = start_sigma;
      planner->method = method;
      planner->T = T;
      planner->noise_level = noise_level;
      planner->initialize();

      while (!planner->finished()) {
        //cout << "Solving optimization" << endl;
        planner->solve(opt_callback);
        //planner->simulate_executions(planner->helper->T);
        planner->simulate_executions(1);
        if (plotting) {
          //emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
          sim_plotter->update_plot_data(&planner->result, &planner->simulated_positions);
        }
      }

      //planner->simulate_executions(planner->helper->T);
      //cout << "final position: " << planner->simulated_positions.back().head<2>().transpose() << endl;
      //cout << "final estimated position: " << planner->start.head<2>().transpose() << endl;
      //cout << "final covariance: " << endl << planner->start_sigma.topLeftCorner<2, 2>() << endl;
      //cout << "goal position: " << planner->goal.head<2>().transpose() << endl;

      cout << endl;

      Vector2d fpos = planner->simulated_positions.back().head<2>();
      Vector2d goalpos = planner->goal.head<2>();
      double dist = (fpos - goalpos).norm();
      cout << "dist: " << dist << endl;
      cout << "goaleps: " << planner->helper->goaleps << endl;
      if (dist < planner->helper->goaleps) {
        cout << "success" << endl;
        ++success;
      }
      cout << "Successful runs: " << success << "/1000" << endl;
    }

    /*
    while (!planner->finished()) {
      //cout << "Solving optimization" << endl;
      planner->solve(opt_callback);
      //planner->simulate_executions(planner->helper->T);
      planner->simulate_executions(1);
      if (plotting) {
        //emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
        sim_plotter->update_plot_data(&planner->result, &planner->simulated_positions);
      }
    }
    */

    cout << "final position: " << planner->simulated_positions.back().head<2>().transpose() << endl;
    cout << "final estimated position: " << planner->start.head<2>().transpose() << endl;
    cout << "final covariance: " << endl << planner->start_sigma.topLeftCorner<2, 2>() << endl;
    cout << "goal position: " << planner->goal.head<2>().transpose() << endl;
    if (plotting) {
      emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
      emit finished_signal();
    }
  }
}

using namespace ToyBSP;

int main(int argc, char *argv[]) {
	bool plotting = true;
	{
		Config config;
		config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
		CommandParser parser(config);
		parser.read(argc, argv, true);
	}
	QApplication app(argc, argv);
	ToyOptimizerTask* task = new ToyOptimizerTask(argc, argv, &app);
	if (plotting) {
		QTimer::singleShot(0, task, SLOT(run_slot()));
		QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
		return app.exec();
	} else {
		task->run();
		return 0;
	}
}
