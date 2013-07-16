#include "bsp/bsp.hpp"
#include "needle.hpp"
#include <deque>
#include <QApplication>
#include <QtCore>
#include <QPolygonF>

using namespace BSP;

namespace NeedleBSP {

  NeedleBSPPlanner::NeedleBSPPlanner() : BSPPlanner<NeedleBSPProblemHelper>() {}

  void NeedleBSPPlanner::initialize() {
    assert(!initialized);
    helper.reset(new NeedleBSPProblemHelper());
    helper->start = start;
    helper->goal = goal;
    helper->start_sigma = start_sigma;
    helper->initialize();
    helper->RRTplan(true);
    helper->T = (int) helper->initial_controls.size();
    state_noise_mean = StateNoiseT::Zero(helper->state_noise_dim);
    state_noise_cov = StateNoiseCovT::Identity(helper->state_noise_dim, helper->state_noise_dim);
    observe_noise_mean = ObserveNoiseT::Zero(helper->observe_noise_dim);
    observe_noise_cov = ObserveNoiseCovT::Identity(helper->observe_noise_dim, helper->observe_noise_dim);
    current_position = start;
    simulated_positions.push_back(current_position);
    initialized = true;
  }

  void NeedleBSPPlanner::custom_simulation_update(StateT* state, VarianceT* sigma, const StateT& actual_state) {
    assert (state != nullptr);
    assert (sigma != nullptr);
    double x = actual_state(0), ultrasound_x = actual_state(3);
    double xmin = ultrasound_x - helper->ultrasound_span * 0.5;
    double xmax = ultrasound_x + helper->ultrasound_span * 0.5;
    Vector2d cur_state = state->head<2>(), new_state;
    Matrix2d cur_sigma = sigma->topLeftCorner<2, 2>(), new_sigma;
    if (x < xmin) {
      truncate_gaussian(Vector2d(1, 0), xmin, cur_state, cur_sigma, &new_state, &new_sigma);
    } else {
      truncate_gaussian(Vector2d(-1, 0), -xmin, cur_state, cur_sigma, &new_state, &new_sigma);
    }
    cout << "state before: " << state->transpose() << endl;
    cout << "sigma before: " << endl << *sigma << endl;
    state->head<2>() = new_state;
    sigma->topLeftCorner<2, 2>() = new_sigma;
    cout << "state: " << state->transpose() << endl;
    cout << "sigma: " << endl << *sigma << endl;
  }

  UltrasoundDeviationCost::UltrasoundDeviationCost(const VarArray& vars, double coeff) : Cost("Ultrasound Deviation"), vars(vars), coeff(coeff) {}

  double UltrasoundDeviationCost::value(const vector<double>& xvec) {
    double val = 0.0;
    for (int i = 0; i < vars.rows(); ++i) {
      StateT state = getVec(xvec, vars.row(i));
      val += fabs(state(0) - state(3));
    }
    return val * coeff;
  }
  
  ConvexObjectivePtr UltrasoundDeviationCost::convex(const vector<double>& xvec, Model* model) {
    ConvexObjectivePtr out(new ConvexObjective(model));
    for (int i = 0; i < vars.rows(); ++i) {
      out->addAbs(exprSub(AffExpr(vars.at(i, 0)), AffExpr(vars.at(i, 3))), coeff);
    }
    return out;
  }


  NeedleBSPProblemHelper::NeedleBSPProblemHelper() : BSPProblemHelper<NeedleBeliefFunc>() {
    input_dt = 0.25;
    needlelen = 0.5;
    goaleps = 0.1;

    ultrasound_span = 0.4;

    set_state_dim(4);
    set_sigma_dof(10);
    set_observe_dim(3);
    set_control_dim(3);
    robot_state_dim = 3;
    robot_control_dim = 2;

    double state_lbs_array[] = {0, 0, -PI*1.25, 0};
    double state_ubs_array[] = {10, 2, PI*1.25, 10};
    double control_lbs_array[] = {-PI*0.25, 0, -3};
    double control_ubs_array[] = {PI*0.25, 2, 3};

    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));

    VarianceT variance_cost = VarianceT::Identity(state_dim, state_dim);
    set_variance_cost(variance_cost);
    set_final_variance_cost(variance_cost * 100);

    ControlCostT control_cost = ControlCostT::Identity(control_dim, control_dim);
    control_cost.bottomRightCorner<1, 1>() = Matrix1d::Identity() * 0.1;//Zero();
    set_control_cost(control_cost * 0.1);//ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void NeedleBSPProblemHelper::configure_problem(OptProb& prob) {
    BSPProblemHelper<NeedleBeliefFunc>::configure_problem(prob);
    add_ultrasound_deviation_cost(prob);
  }

  void NeedleBSPProblemHelper::add_ultrasound_deviation_cost(OptProb& prob) {
    prob.addCost(CostPtr(new UltrasoundDeviationCost(state_vars, 0.1)));
  }

  void NeedleBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      // box constraint of buffer goaleps around goal position
      Var optvar = state_vars.at(T, i);
      prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i)+goaleps) ), INEQ);
      prob.addLinearConstraint(exprSub(exprSub(AffExpr(0), AffExpr(optvar)), (-goal(i)+goaleps) ), INEQ);
    }
  }

  void NeedleBSPProblemHelper::RRTplan(bool compute) {
    if (compute) {

      srand(time(0));

      vector<RRTNode> rrtTree;
      RRTNode startNode;
      startNode.x = start.head<3>();
      rrtTree.push_back(startNode);

      Vector2d poserr = (startNode.x.head<2>() - goal.head<2>());

      double state_lbs_array[] = {0, 0, -PI*1.25};
      double state_ubs_array[] = {10, 2, PI*1.25};
      double control_lbs_array[] = {-PI*0.25, 0};
      double control_ubs_array[] = {PI*0.25, 3};

      int numiter = 0;
      while (poserr.squaredNorm() > goaleps*goaleps || numiter < 100) {

        cout << ".";
        RobotStateT sample;
        for (int xd = 0; xd < robot_state_dim; ++xd) {
          sample(xd) = (((double) rand()) / RAND_MAX) * (state_ubs_array[xd] - state_lbs_array[xd]) + state_lbs_array[xd];
        }

        // Check sample for collisions, turned off for now

        int node = -1;
        double mindist = 9e99;
        for (int j = 0; j < (int) rrtTree.size(); ++j) {
          double ddist = (rrtTree[j].x - sample).squaredNorm();
          if (ddist < mindist) {
            mindist = ddist;
            node = j;
          }
        }
        if (node == -1) {
          continue;
        }

        RobotControlT input;
        for (int ud = 0; ud < robot_control_dim; ++ud) {
          input(ud) = (((double) rand()) / RAND_MAX) * (control_ubs_array[ud] - control_lbs_array[ud]) + control_lbs_array[ud];
        }

        StateNoiseT zero_state_noise = StateNoiseT::Zero(state_noise_dim);
        RobotStateT new_x = state_func->call((StateT) concat(rrtTree[node].x, Vector4d::Zero()), (ControlT) concat(input, Vector4d::Zero()), zero_state_noise).head<3>();

        bool valid = true;
        for (int xd = 0; xd < robot_state_dim; ++xd) {
          if (new_x(xd) < state_lbs[xd] || new_x(xd) > state_ubs[xd]) {
            valid = false;
            break;
          }
        }
        if (!valid) {
          continue;
        }

        // Check edge for collisions here, turned off for now

        RRTNode newnode;
        newnode.x = new_x;
        newnode.u = input;
        newnode.bp = node;

        rrtTree.push_back(newnode);
        Vector4d edge;
        edge << rrtTree[node].x(0), rrtTree[node].x(1), newnode.x(0), newnode.x(1);
        rrt_edges.push_back(edge);

        poserr = (newnode.x.head<2>() - goal.head<2>());
        ++numiter;
      }
      cout << endl;

      deque<RRTNode> path;

      int i = (int)rrtTree.size() - 1;
      RRTNode node;
      node.x = rrtTree[i].x;
      node.u = RobotControlT::Zero(robot_control_dim);
      path.push_front(node);

      while (i != 0) {
        node.u = rrtTree[i].u;
        i = rrtTree[i].bp;
        node.x = rrtTree[i].x;
        node.bp = -1;

        path.push_front(node);
      }

      initial_controls.clear();
      for (int i = 0; i < (int)path.size()-1; ++i) {
        initial_controls.push_back((ControlT) concat(path[i].u, Vector4d::Zero()));
        cout << path[i].u(0) << " " << path[i].u(1) << endl;
      }
      cout << "T: " << initial_controls.size() << endl;
    } else {

      ifstream fptr("needle-rrt-seq.txt", ios::in);
      if (!fptr.is_open()) {
        cerr << "Could not open file, check location" << endl;
        std::exit(-1);
      }
      int nu;
      fptr >> nu;
      //cout << "nu: " << nu << endl;
      initial_controls.clear();
      ControlT u = ControlT::Zero(control_dim);
      for (int i = 0; i < nu; ++i) {
        fptr >> u(0) >> u(1);
        initial_controls.push_back(u);
        //cout << u(0) << " " << u(1) << endl;
      }
      if (fptr.is_open()) fptr.close();
    }

  }

  NeedleStateFunc::NeedleStateFunc() : StateFunc<StateT, ControlT, StateNoiseT>() {}

  NeedleStateFunc::NeedleStateFunc(BSPProblemHelperBasePtr helper) :
                  StateFunc<StateT, ControlT, StateNoiseT>(helper), needle_helper(boost::static_pointer_cast<NeedleBSPProblemHelper>(helper)) {}

  StateT NeedleStateFunc::operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const {
    StateT new_x(state_dim);
    /* Euler integration */
    //new_x(0) = x(0) + needle_helper->input_dt * x(3) * cos(x(2));
    //new_x(1) = x(1) + needle_helper->input_dt * x(3) * sin(x(2));
    //new_x(2) = x(2) + needle_helper->input_dt * x(3) * tan(u(0)) / needle_helper->needlelen;
    //new_x(3) = x(3) + needle_helper->input_dt * u(1);
    /* RK4 integration */
    StateT xtmp(state_dim), x1(state_dim), x2(state_dim), x3(state_dim), x4(state_dim);
    xtmp = x;
    x1(0) = u(1) * cos(xtmp(2));
    x1(1) = u(1) * sin(xtmp(2));
    x1(2) = u(1) * tan(u(0)) / needle_helper->needlelen;
    xtmp = x + 0.5*needle_helper->input_dt*x1;
    x2(0) = u(1) * cos(xtmp(2));
    x2(1) = u(1) * sin(xtmp(2));
    x2(2) = u(1) * tan(u(0)) / needle_helper->needlelen;
    xtmp = x + 0.5*needle_helper->input_dt*x2;
    x3(0) = u(1) * cos(xtmp(2));
    x3(1) = u(1) * sin(xtmp(2));
    x3(2) = u(1) * tan(u(0)) / needle_helper->needlelen;
    xtmp = x + needle_helper->input_dt*x3;
    x4(0) = u(1) * cos(xtmp(2));
    x4(1) = u(1) * sin(xtmp(2));
    x4(2) = u(1) * tan(u(0)) / needle_helper->needlelen;

    new_x = x + needle_helper->input_dt/6.0*(x1+2.0*(x2+x3)+x4);
    new_x(3) = x(3) + u(2) * needle_helper->input_dt;
    new_x.head<3>() += 0.1 * m.head<3>();
    new_x(3) += 0.01 * m(3);

    double umag = u.squaredNorm();
    //return new_x + 0.1*umag*m;
    return new_x;// + 0.01*m;
  }

  NeedleObserveFunc::NeedleObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  NeedleObserveFunc::NeedleObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), needle_helper(boost::static_pointer_cast<NeedleBSPProblemHelper>(helper)) {}

  ObserveT NeedleObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT ret(observe_dim);
    ret(0) = x(0);
    ret(1) = x(1);
    ret(2) = x(3);
    return ret + 0.01 * n;
  }

  bool NeedleObserveFunc::sgndist(double x, double ultrasound_x, Vector1d* dists) const {
    double xmin = ultrasound_x - needle_helper->ultrasound_span * 0.5;
    double xmax = ultrasound_x + needle_helper->ultrasound_span * 0.5;
    if (xmin <= x && x <= xmax) {
      (*dists)(0) = -fmin(x - xmin, xmax - x);
    } else {
      (*dists)(0) = fmin(fabs(x - xmin), fabs(xmax - x));
    }
    return (*dists)(0) < 0;
  }

  ObserveT NeedleObserveFunc::observation_masks(const StateT& x, double approx_factor) const {
    ObserveT ret(observe_dim);
    Vector1d dists;
    sgndist(x(0), x(3), &dists);
    double minval = 1e-4;
    double tol = 0.1;

    ret(2) = 1;
    if (approx_factor < 0) {
      ret(0) = ret(1) = dists(0) <= 0 ? 1 : 0;
    } else {
      ret(0) = ret(1) = 1. - sigmoid(approx_factor * (dists(0) + tol));
    }
  }

  NeedleBeliefFunc::NeedleBeliefFunc() : EkfBeliefFunc<NeedleStateFunc, NeedleObserveFunc, BeliefT>() {
    this->approx_factor = 1;
  }

  NeedleBeliefFunc::NeedleBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
   EkfBeliefFunc<NeedleStateFunc, NeedleObserveFunc, BeliefT>(helper, f, h), needle_helper(boost::static_pointer_cast<NeedleBSPProblemHelper>(helper)) {
    this->approx_factor = 1;
  }

  NeedlePlotter::NeedlePlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   needle_helper(boost::static_pointer_cast<NeedleBSPProblemHelper>(helper)) {}

  void NeedlePlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
    assert(distmap != nullptr);
    for (int j = 0; j < distmap->height(); ++j) {
      QRgb *line = (QRgb*) distmap->scanLine(j);
      for (int i = 0; i < distmap->width(); ++i) {
        if (state == nullptr) {
          line[i] = qRgb(0, 0, 0);
        } else {
          double x = unscale_x(i);
          Vector1d dists;
          needle_helper->belief_func->h->sgndist(x, (*state)(3), &dists);
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

  NeedleOptPlotter::NeedleOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   NeedlePlotter(x_min, x_max, y_min, y_max, helper, parent),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void NeedleOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
      // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      if (states_opt.size() > 0) {
        compute_distmap(&distmap, &states_opt[0], needle_helper->belief_func->approx_factor);
      } else {
        compute_distmap(&distmap, nullptr, needle_helper->belief_func->approx_factor);
      }
    }
    painter.drawImage(0, 0, distmap);
    QPen cvx_cov_pen(Qt::red, 3, Qt::SolidLine);
    QPen path_pen(Qt::red, 3, Qt::SolidLine);
    QPen pos_pen(Qt::red, 8, Qt::SolidLine);
    QPen ultrasound_path_pen(Qt::white, 3, Qt::SolidLine);
    QPen ultrasound_pos_pen(Qt::white, 8, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(cvx_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].head<2>(), sigmas_opt[i].topLeftCorner(2,2), painter, 0.5);
    }

    painter.setPen(path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](0), states_opt[i](1), states_opt[i+1](0), states_opt[i+1](1), painter);
    }
    //vector<Vector4d>& edges = needle_helper->rrt_edges;
    //for (int i = 0; i < edges.size(); ++i) {
    //  draw_line(edges[i](0), edges[i](1), edges[i](2), edges[i](3), painter);
    //}
    painter.setPen(pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](0), states_opt[i](1), painter);
    }

    painter.setPen(ultrasound_path_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](3), -1, painter);
    }
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](3), -1, states_opt[i+1](3), -1, painter);
    }


    // draw beliefs computed using belief dynamics
    QPen cvx_cov_pen2(Qt::blue, 3, Qt::SolidLine);
    painter.setPen(cvx_cov_pen2);
    for (int i = 0; i < states_actual.size(); ++i) {
      draw_ellipse(states_actual[i].head<2>(), sigmas_actual[i].topLeftCorner(2,2), painter, 0.5);
    }
    QPen pos_pen2(Qt::blue, 8, Qt::SolidLine);
    painter.setPen(pos_pen2);
    for (int i = 0; i < states_actual.size(); ++i) {
      draw_point(states_actual[i](0), states_actual[i](1), painter);
    }
  }

  void NeedleOptPlotter::update_plot_data(void* data) {
    vector<StateT> new_states_opt, new_states_actual;
    vector<VarianceT> new_sigmas_opt, new_sigmas_actual;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = needle_helper->belief_func->approx_factor;
    DblVec* xvec = (DblVec* ) data;
    BeliefT cur_belief_opt, cur_belief_actual;
    needle_helper->belief_func->compose_belief(needle_helper->start, matrix_sqrt(needle_helper->start_sigma), &cur_belief_actual);
    int T = needle_helper->get_T();

    for (int i = 0; i <= T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;

      cur_belief_opt = (BeliefT) getVec(*xvec, needle_helper->belief_vars.row(i));
      needle_helper->belief_func->extract_state(cur_belief_opt, &cur_state);
      needle_helper->belief_func->extract_sigma(cur_belief_opt, &cur_sigma);
      new_states_opt.push_back(cur_state);
      new_sigmas_opt.push_back(cur_sigma);

      needle_helper->belief_func->extract_state(cur_belief_actual, &cur_state);
      needle_helper->belief_func->extract_sigma(cur_belief_actual, &cur_sigma);
      new_states_actual.push_back(cur_state);
      new_sigmas_actual.push_back(cur_sigma);
      if (i < T) {
        cur_belief_actual = needle_helper->belief_func->call(cur_belief_actual, (ControlT) getVec(*xvec, needle_helper->control_vars.row(i)));
      }
    }
    states_opt = new_states_opt; states_actual = new_states_actual;
    sigmas_opt = new_sigmas_opt; sigmas_actual = new_sigmas_actual;

    needle_helper->belief_func->approx_factor = cur_approx_factor;
    this->repaint();
  }

  NeedleSimulationPlotter::NeedleSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   NeedlePlotter(x_min, x_max, y_min, y_max, helper, parent),
   distmap(400, 400, QImage::Format_RGB32) {}

  void NeedleSimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    //distmap = QImage(width(), height(), QImage::Format_RGB32);
    //if (simulated_positions.size() <= 0) {
    //  return;
    //}
    //compute_distmap(&distmap, &simulated_positions.back(), -1);
    //painter.drawImage(0, 0, distmap);
    QPen sim_prev_pos_pen(QColor(255, 0, 0, 100), 4, Qt::SolidLine);
    QPen sim_cur_pos_pen(QColor(255, 0, 0, 255), 4, Qt::SolidLine);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(sim_cur_pos_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i+1](0), simulated_positions[i+1](1), painter);
    }

  }

  void NeedleSimulationPlotter::update_plot_data(void* data_x, void* data_sim) {
    simulated_positions = *((vector<StateT>*) data_sim);
    DblVec* xvec = (DblVec* ) data_x;
    vector<StateT> new_states;
    vector<VarianceT> new_sigmas;
    BeliefT cur_belief;
    needle_helper->belief_func->compose_belief(needle_helper->start, matrix_sqrt(needle_helper->start_sigma), &cur_belief);
    for (int i = 0; i <= needle_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      needle_helper->belief_func->extract_state(cur_belief, &cur_state);
      needle_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < needle_helper->T) cur_belief = needle_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, needle_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;

    this->repaint();
  }

  NeedleOptimizerTask::NeedleOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  NeedleOptimizerTask::NeedleOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void NeedleOptimizerTask::stage_plot_callback(boost::shared_ptr<NeedleOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
    //wait_to_proceed(boost::bind(&NeedleOptPlotter::update_plot_data, plotter, &x));
  }


  void NeedleOptimizerTask::run() {
    bool plotting = true;

    double start_vec_array[] = {0, 1, 0, 0};
    double goal_vec_array[] = {10, 1.5, 0, 0};

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
    Matrix4d start_sigma = Matrix4d::Identity()*1;
    //start_sigma.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 2;

    NeedleBSPPlannerPtr planner(new NeedleBSPPlanner());
    planner->start = start;
    planner->goal = goal;
    planner->start_sigma = start_sigma;
    planner->initialize();

    boost::shared_ptr<NeedleSimulationPlotter> sim_plotter;
    boost::shared_ptr<NeedleOptPlotter> opt_plotter;

    if (plotting) {
      double x_min = -2, x_max = 12, y_min = -2, y_max = 4;
      sim_plotter.reset(create_plotter<NeedleSimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      opt_plotter.reset(create_plotter<NeedleOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      opt_plotter->show();
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    if (plotting) {
      opt_callback = boost::bind(&NeedleOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    while (!planner->finished()) {
      planner->solve(opt_callback);
      planner->simulate_execution();//s(planner->helper->T);
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
      }
    }

    emit finished_signal();
    
  }
}

using namespace NeedleBSP;

int main(int argc, char *argv[]) {
	QApplication app(argc, argv);
	NeedleOptimizerTask* task = new NeedleOptimizerTask(argc, argv, &app);
	QTimer::singleShot(0, task, SLOT(run_slot()));
	QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
	return app.exec();
}
