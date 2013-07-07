#include "bsp/bsp.hpp"
#include "car.hpp"
#include <deque>
#include <QApplication>
#include <QtCore>
#include <QPolygonF>

using namespace BSP;

namespace CarBSP {

  CarBSPPlanner::CarBSPPlanner() : BSPPlanner<CarBSPProblemHelper>() {}

  void CarBSPPlanner::initialize() {
    assert(!initialized);
    helper.reset(new CarBSPProblemHelper());
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

  CarBSPProblemHelper::CarBSPProblemHelper() : BSPProblemHelper<CarBeliefFunc>() {
    input_dt = 0.25;
    carlen = 0.5;
    goaleps = 0.1;

    car_camera_depth = 2;
    car_camera_span_angle = PI / 3;

    set_state_dim(7);
    set_sigma_dof(28);
    set_observe_dim(4);
    set_control_dim(6);
    robot_state_dim = 3;
    robot_control_dim = 2;

    double state_lbs_array[] = {-10, -10, -PI*1.25, -10, -10, -10, -10};
    double state_ubs_array[] = {10, 10, PI*1.25, 10, 10, 10, 10};
    double control_lbs_array[] = {-PI*0.25, 0, -10, -10, -10, -10};
    double control_ubs_array[] = {PI*0.25, 3, 10, 10, 10, 10};

    // static controls:
    //double control_lbs_array[] = {-PI*0.25, 0, 0, 0, 0, 0};
    //double control_ubs_array[] = {PI*0.25, 3, 0, 0, 0, 0};

    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));

    VarianceT variance_cost = VarianceT::Identity(state_dim, state_dim);
    //variance_cost.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 0.01;
    set_variance_cost(variance_cost);//VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(variance_cost * 100);//VarianceT::Identity(state_dim, state_dim)*100);

    ControlCostT control_cost = ControlCostT::Identity(control_dim, control_dim);
    //control_cost.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 0.1;//Zero();
    set_control_cost(control_cost * 0.1);//ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void CarBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      // box constraint of buffer goaleps around goal position
      Var optvar = state_vars.at(T, i);
      prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i)+goaleps) ), INEQ);
      prob.addLinearConstraint(exprSub(exprSub(AffExpr(0), AffExpr(optvar)), (-goal(i)+goaleps) ), INEQ);
    }
  }

  void CarBSPProblemHelper::RRTplan(bool compute) {
    if (compute) {

      srand(time(0));

      vector<RRTNode> rrtTree;
      RRTNode startNode;
      startNode.x = start.head<3>();
      rrtTree.push_back(startNode);

      Vector2d poserr = (startNode.x.head<2>() - goal.head<2>());

      double state_lbs_array[] = {-6, -6, -PI};
      double state_ubs_array[] = {6, 6, PI};
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

      ifstream fptr("car-rrt-seq.txt", ios::in);
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
    x1(0) = u(1) * cos(xtmp(2));
    x1(1) = u(1) * sin(xtmp(2));
    x1(2) = u(1) * tan(u(0)) / car_helper->carlen;
    xtmp = x + 0.5*car_helper->input_dt*x1;
    x2(0) = u(1) * cos(xtmp(2));
    x2(1) = u(1) * sin(xtmp(2));
    x2(2) = u(1) * tan(u(0)) / car_helper->carlen;
    xtmp = x + 0.5*car_helper->input_dt*x2;
    x3(0) = u(1) * cos(xtmp(2));
    x3(1) = u(1) * sin(xtmp(2));
    x3(2) = u(1) * tan(u(0)) / car_helper->carlen;
    xtmp = x + car_helper->input_dt*x3;
    x4(0) = u(1) * cos(xtmp(2));
    x4(1) = u(1) * sin(xtmp(2));
    x4(2) = u(1) * tan(u(0)) / car_helper->carlen;

    new_x = x + car_helper->input_dt/6.0*(x1+2.0*(x2+x3)+x4);
    new_x.tail<4>() = x.tail<4>() + u.tail<4>() * car_helper->input_dt;
    new_x.head<3>() += 0.01 * m.head<3>();
    new_x.tail<4>() += 0.01 * m.tail<4>();

    double umag = u.squaredNorm();
    //return new_x + 0.1*umag*m;
    return new_x;// + 0.01*m;
  }

  CarObserveFunc::CarObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  CarObserveFunc::CarObserveFunc(BSPProblemHelperBasePtr helper) :
    ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  ObserveT CarObserveFunc::unnoisy_observation(const StateT& x) const {
    ObserveT ret(observe_dim);
    Vector2d car_pos = x.head<2>(), l1 = x.middleRows<2>(3), l2 = x.middleRows<2>(5);
    double car_angle = x(2);
    //ret(0) = car_velocity;
    //ret(1) = car_pos.x() - l1.x();
    ret(0) = (car_pos - l1).norm();
    ret(1) = atan2(car_pos.y() - l1.y(), car_pos.x() - l1.x()) - car_angle;
    ret(2) = (car_pos - l2).norm();
    ret(3) = atan2(car_pos.y() - l2.y(), car_pos.x() - l2.x()) - car_angle;
    return ret;
  }

  ObserveT CarObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    return unnoisy_observation(x) + 0.01 * compute_inverse_gamma(x, -1) * n;
  }

  ObserveT CarObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n, double approx_factor) const {
    return unnoisy_observation(x) + 0.01 * compute_inverse_gamma(x, approx_factor) * n;
  }

  Beam2D robot_fov(double x, double y, double angle, double camera_depth, double camera_span_angle) {
    Beam2D out;
    double x0 = x, y0 = y;
    double skew_depth = camera_depth / cos(camera_span_angle * 0.5);
    double x1 = x0 + skew_depth * cos(angle - camera_span_angle * 0.5),
           y1 = y0 + skew_depth * sin(angle - camera_span_angle * 0.5);
    double x2 = x0 + skew_depth * cos(angle + camera_span_angle * 0.5),
           y2 = y0 + skew_depth * sin(angle + camera_span_angle * 0.5);
    out.base = Point(x0, y0);
    out.a = Point(x1, y1);
    out.b = Point(x2, y2);
    return out;
  }

  bool CarObserveFunc::sgndist(const StateT& x, Vector2d* dists) const {
    Beam2D fov = robot_fov(x(0), x(1), x(2), car_helper->car_camera_depth, car_helper->car_camera_span_angle);
    //cout << "fov: " << fov << endl;
    (*dists)(0) = Geometry2D::sgndist(x.middleRows<2>(3), fov);
    (*dists)(1) = Geometry2D::sgndist(x.middleRows<2>(5), fov);
    //cout << "dists: " << dists->transpose() << endl;
    return (*dists)(0) < 0 || (*dists)(1) < 0;
  }

  ObserveMatT CarObserveFunc::compute_gamma(const StateT& x, double approx_factor) const {
    Vector2d dists;
    sgndist(x, &dists);
    double tol = 0.1;
    double gamma1, gamma2;
    if (approx_factor < 0) {
      gamma1 = dists(0) <= 0 ? 1 : 0; 
      gamma2 = dists(1) <= 0 ? 1 : 0;
    } else {
      gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
    }
    ObserveMatT gamma(observe_dim, observe_dim);

    gamma << gamma1,      0,      0,      0,
                  0, gamma1,      0,      0,
                  0,      0, gamma2,      0,
                  0,      0,      0, gamma2;

    return gamma;
  }

  ObserveMatT CarObserveFunc::compute_inverse_gamma(const StateT& x, double approx_factor) const {
    Vector2d dists;
    sgndist(x, &dists);
    double minval = 1e-4;
    double tol = 0.1;
    double invgamma1, invgamma2;
    if (approx_factor < 0) {
      invgamma1 = dists(0) <= 0 ? 1 : 1/minval;
      invgamma2 = dists(1) <= 0 ? 1 : 1/minval;
    } else {
      double gamma1 = 1. - (1./(1.+exp(-approx_factor*(dists(0)+tol))));
      double gamma2 = 1. - (1./(1.+exp(-approx_factor*(dists(1)+tol))));
      invgamma1 = 1. / fmax(gamma1, minval);
      invgamma2 = 1. / fmax(gamma2, minval);
    }

    ObserveMatT invgamma(observe_dim, observe_dim);
    invgamma << invgamma1,         0,         0,         0,
                        0, invgamma1,         0,         0,
                        0,         0, invgamma2,         0,
                        0,         0,         0, invgamma2;

    return invgamma;
  }

  CarBeliefFunc::CarBeliefFunc() : BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT>() {
    this->approx_factor = 1;
  }

  CarBeliefFunc::CarBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h) :
   BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT>(helper, f, h), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {
    this->approx_factor = 1;
  }

  CarPlotter::CarPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   BSPQtPlotter(x_min, x_max, y_min, y_max, helper, parent),
   ProblemState(helper),
   car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  void CarPlotter::compute_distmap(QImage* distmap, StateT* state, double approx_factor) {
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
          Vector2d dists;
          car_helper->belief_func->h->sgndist((StateT) concat(pos, state->tail<5>()), &dists);
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
  }

  void CarPlotter::draw_beam_2d(const Beam2D& beam, QPainter& painter) {
    QPolygonF polygon;
    polygon << QPointF(scale_x(beam.base.x()), scale_y(beam.base.y()));
    polygon << QPointF(scale_x(beam.a.x()), scale_y(beam.a.y()));
    polygon << QPointF(scale_x(beam.b.x()), scale_y(beam.b.y()));
    painter.drawPolygon(polygon);
  }

  CarOptPlotter::CarOptPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   CarPlotter(x_min, x_max, y_min, y_max, helper, parent),
   old_approx_factor(-1), cur_approx_factor(-1), distmap(400, 400, QImage::Format_RGB32) {}

  void CarOptPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    //if (cur_approx_factor != old_approx_factor || distmap.height() != height() || distmap.width() != width()) {
    //  // replot distmap
    //  distmap = QImage(width(), height(), QImage::Format_RGB32);
    //  if (states_opt.size() > 0) {
    //    compute_distmap(&distmap, &states_opt[0], car_helper->belief_func->approx_factor);
    //  } else {
    //    compute_distmap(&distmap, NULL, car_helper->belief_func->approx_factor);
    //  }
    //}
    //painter.drawImage(0, 0, distmap);
    QPen cvx_cov_pen(Qt::red, 3, Qt::SolidLine);
    QPen path_pen(Qt::red, 3, Qt::SolidLine);
    QPen pos_pen(Qt::red, 8, Qt::SolidLine);
    QPen fov_pen(QColor(255, 255, 255, 30), 1, Qt::SolidLine);
    QBrush fov_brush(QColor(255, 255, 255, 30));
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(cvx_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].head<2>(), sigmas_opt[i].topLeftCorner(2,2), painter, 0.5);
    }
    painter.setPen(fov_pen);
    QBrush prev_brush = painter.brush();
    painter.setBrush(fov_brush);

    for (int i = 0; i < states_opt.size(); ++i) {
      Beam2D fov = robot_fov(states_opt[i](0), states_opt[i](1), states_opt[i](2), car_helper->car_camera_depth, car_helper->car_camera_span_angle);
      cout << "fov: " << fov << endl;
      draw_beam_2d(fov, painter);
    }
    painter.setBrush(prev_brush);

    painter.setPen(path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](0), states_opt[i](1), states_opt[i+1](0), states_opt[i+1](1), painter);
    }
    //vector<Vector4d>& edges = car_helper->rrt_edges;
    //for (int i = 0; i < edges.size(); ++i) {
    //  draw_line(edges[i](0), edges[i](1), edges[i](2), edges[i](3), painter);
    //}
    painter.setPen(pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](0), states_opt[i](1), painter);
    }

    QPen sensor_cov_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].middleRows<2>(3), sigmas_opt[i].block<2, 2>(3, 3), painter, 0.5);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](3), states_opt[i](4), states_opt[i+1](3), states_opt[i+1](4), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](3), states_opt[i](4), painter);
    }

    painter.setPen(sensor_cov_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_opt[i].middleRows<2>(5), sigmas_opt[i].block<2, 2>(5, 5), painter, 0.5);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < states_opt.size() - 1; ++i) {
      draw_line(states_opt[i](5), states_opt[i](6), states_opt[i+1](5), states_opt[i+1](6), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_opt.size(); ++i) {
      draw_point(states_opt[i](5), states_opt[i](6), painter);
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

  void CarOptPlotter::update_plot_data(void* data) {
    vector<StateT> new_states_opt, new_states_actual;
    vector<VarianceT> new_sigmas_opt, new_sigmas_actual;
    old_approx_factor = cur_approx_factor;
    cur_approx_factor = car_helper->belief_func->approx_factor;
    DblVec* xvec = (DblVec* ) data;
    BeliefT cur_belief_opt, cur_belief_actual;
    car_helper->belief_func->compose_belief(car_helper->start, matrix_sqrt(car_helper->start_sigma), &cur_belief_actual);
    int T = car_helper->get_T();

    for (int i = 0; i <= T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;

      cur_belief_opt = (BeliefT) getVec(*xvec, car_helper->belief_vars.row(i));
      car_helper->belief_func->extract_state(cur_belief_opt, &cur_state);
      car_helper->belief_func->extract_sigma(cur_belief_opt, &cur_sigma);
      new_states_opt.push_back(cur_state);
      new_sigmas_opt.push_back(cur_sigma);

      car_helper->belief_func->extract_state(cur_belief_actual, &cur_state);
      car_helper->belief_func->extract_sigma(cur_belief_actual, &cur_sigma);
      new_states_actual.push_back(cur_state);
      new_sigmas_actual.push_back(cur_sigma);
      if (i < T) {
        cur_belief_actual = car_helper->belief_func->call(cur_belief_actual, (ControlT) getVec(*xvec, car_helper->control_vars.row(i)));
      }
    }
    states_opt = new_states_opt; states_actual = new_states_actual;
    sigmas_opt = new_sigmas_opt; sigmas_actual = new_sigmas_actual;

    car_helper->belief_func->approx_factor = cur_approx_factor;
    this->repaint();
  }

  CarSimulationPlotter::CarSimulationPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent) :
   CarPlotter(x_min, x_max, y_min, y_max, helper, parent),
   distmap(400, 400, QImage::Format_RGB32) {}

  void CarSimulationPlotter::paintEvent(QPaintEvent* ) {
    QPainter painter(this);
    //distmap = QImage(width(), height(), QImage::Format_RGB32);
    //if (simulated_positions.size() <= 0) {
    //  return;
    //}
    //compute_distmap(&distmap, &simulated_positions.back(), -1);
    //painter.drawImage(0, 0, distmap);
    QPen sim_prev_pos_pen(QColor(255, 0, 0, 100), 4, Qt::SolidLine);
    QPen sim_cur_pos_pen(QColor(255, 0, 0, 255), 4, Qt::SolidLine);
    QPen fov_pen(QColor(255, 255, 255, 30), 1, Qt::SolidLine);
    QBrush fov_brush(QColor(255, 255, 255, 30));
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    painter.setPen(sim_cur_pos_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i+1](0), simulated_positions[i+1](1), painter);
    }

    QBrush prev_brush = painter.brush();
    painter.setPen(fov_pen);
    painter.setBrush(fov_brush);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      Beam2D fov = robot_fov(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i](2), car_helper->car_camera_depth, car_helper->car_camera_span_angle);
      draw_beam_2d(fov, painter);
    }
    painter.setBrush(prev_brush);

    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](3), simulated_positions[i](4), simulated_positions[i+1](3), simulated_positions[i+1](4), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](3), simulated_positions[i](4), painter);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](5), simulated_positions[i](6), simulated_positions[i+1](5), simulated_positions[i+1](6), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](5), simulated_positions[i](6), painter);
    }
  }

  void CarSimulationPlotter::update_plot_data(void* data_x, void* data_sim) {
    simulated_positions = *((vector<StateT>*) data_sim);
    DblVec* xvec = (DblVec* ) data_x;
    vector<StateT> new_states;
    vector<VarianceT> new_sigmas;
    BeliefT cur_belief;
    car_helper->belief_func->compose_belief(car_helper->start, matrix_sqrt(car_helper->start_sigma), &cur_belief);
    for (int i = 0; i <= car_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      car_helper->belief_func->extract_state(cur_belief, &cur_state);
      car_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < car_helper->T) cur_belief = car_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, car_helper->control_vars.row(i)));
    }
    states = new_states;
    sigmas = new_sigmas;

    this->repaint();
  }

  CarOptimizerTask::CarOptimizerTask(QObject* parent) : BSPOptimizerTask(parent) {}

  CarOptimizerTask::CarOptimizerTask(int argc, char **argv, QObject* parent) : BSPOptimizerTask(argc, argv, parent) {}

  void CarOptimizerTask::stage_plot_callback(boost::shared_ptr<CarOptPlotter> plotter, OptProb*, DblVec& x) {
    plotter->update_plot_data(&x);
    //wait_to_proceed(boost::bind(&CarOptPlotter::update_plot_data, plotter, &x));
  }


  void CarOptimizerTask::run() {
    bool plotting = true;

    double start_vec_array[] = {-5, 2, -PI*0.5, 0, 2, 0, -2};
    double goal_vec_array[] = {-5, -2, 0, 0, 0, 0, 0};

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

    Vector7d start = toVectorXd(start_vec);
    Vector7d goal = toVectorXd(goal_vec);
    Matrix7d start_sigma = Matrix7d::Identity()*1;
    start_sigma.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 2;

    CarBSPPlannerPtr planner(new CarBSPPlanner());
    planner->start = start;
    planner->goal = goal;
    planner->start_sigma = start_sigma;
    planner->initialize();

    boost::shared_ptr<CarSimulationPlotter> sim_plotter;
    boost::shared_ptr<CarOptPlotter> opt_plotter;

    if (plotting) {
      double x_min = -7, x_max = 7, y_min = -7, y_max = 7;
      sim_plotter.reset(create_plotter<CarSimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      opt_plotter.reset(create_plotter<CarOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      opt_plotter->show();
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    if (plotting) {
      opt_callback = boost::bind(&CarOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    while (!planner->finished()) {
      planner->solve(opt_callback);
      planner->simulate_executions(planner->helper->T);
      if (plotting) {
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
      }
    }

    emit finished_signal();
    
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
