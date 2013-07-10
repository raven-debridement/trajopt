#include "bsp/bsp.hpp"
#include "car_static.hpp"
#include <deque>
#include <QApplication>
#include <QtCore>
#include <QPolygonF>


using namespace BSP;

namespace CarBSP {

  const static char car_rrt_filename[] = "car-rrt-seq.txt";

  CarBSPPlanner::CarBSPPlanner() : BSPPlanner<CarBSPProblemHelper>() {}

  void CarBSPPlanner::initialize() {
    BSPPlanner<CarBSPProblemHelper>::initialize();
  }

  void CarBSPPlanner::initialize_optimizer_parameters(BSPTrustRegionSQP& opt) {
    opt.max_iter_                   = 250;
    opt.merit_error_coeff_          = 100;
    opt.merit_coeff_increase_ratio_ = 10;
    opt.max_merit_coeff_increases_  = 2;
    opt.trust_shrink_ratio_         = 0.8;
    opt.trust_expand_ratio_         = 1.2;
    opt.min_trust_box_size_         = 0.001;
    opt.min_approx_improve_         = 1e-2;
    opt.min_approx_improve_frac_    = 1e-4;
    opt.improve_ratio_threshold_    = 0.1;
    opt.trust_box_size_             = 1;
    opt.cnt_tolerance_              = 1e-06;
  }

  void CarBSPPlanner::custom_simulation_update(StateT* state, VarianceT* sigma) {
    assert (state != NULL);
    assert (sigma != NULL);
    if ((*state)(2) > 2*PI) {
      (*state)(2) = 2*PI;
    } else if ((*state)(2) < -2*PI) {
      (*state)(2) = -2*PI;
    }
  }


  CarBSPProblemHelper::CarBSPProblemHelper() : BSPProblemHelper<CarBeliefFunc>() {
    input_dt = 0.25;
    carlen = 0.5;
    goaleps = 0.1;

    car_camera_depth = 2.5;
    car_camera_span_angle = PI / 4;

    set_state_dim(7);
    set_sigma_dof(28);
    set_observe_dim(4);
    set_control_dim(2);
    robot_state_dim = 3;
    robot_control_dim = 2;

    double state_lbs_array[] = {-10, -10, -PI*2, -10, -10, -10, -10};
    double state_ubs_array[] = {10, 10, PI*2, 10, 10, 10, 10};

    // static controls:
    double control_lbs_array[] = {-PI*0.25, 0};//, 0, 0, 0, 0};
    double control_ubs_array[] = {PI*0.25, 5};//, 0, 0, 0, 0};//0.00001, 0.00001, 0.00001, 0.00001};

    set_state_bounds(DblVec(state_lbs_array, end(state_lbs_array)), DblVec(state_ubs_array, end(state_ubs_array)));
    set_control_bounds(DblVec(control_lbs_array, end(control_lbs_array)), DblVec(control_ubs_array, end(control_ubs_array)));

    VarianceT variance_cost = VarianceT::Identity(state_dim, state_dim);
    variance_cost.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 0;
    set_variance_cost(variance_cost);//VarianceT::Identity(state_dim, state_dim));
    set_final_variance_cost(variance_cost * 100);//VarianceT::Identity(state_dim, state_dim)*100);

    ControlCostT control_cost = ControlCostT::Identity(control_dim, control_dim);
    //control_cost.bottomRightCorner<4, 4>() = Matrix4d::Identity() * 0;//.1;//Zero();
    set_control_cost(control_cost * 0);//0.01);//ControlCostT::Identity(control_dim, control_dim)*0.1);
  }

  void CarBSPProblemHelper::add_goal_constraint(OptProb& prob) {
    for (int i = 0; i < 2; ++i) {
      // box constraint of buffer goaleps around goal position
      Var optvar = state_vars.at(T, i);
      //prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i))), INEQ);
      prob.addLinearConstraint(exprSub(AffExpr(optvar), (goal(i)+goaleps) ), INEQ);
      prob.addLinearConstraint(exprSub(exprSub(AffExpr(0), AffExpr(optvar)), (-goal(i)+goaleps) ), INEQ);
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
    x1(0) = u(1) * cos((double)xtmp(2));
    x1(1) = u(1) * sin((double)xtmp(2));
    x1(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
    xtmp = x + 0.5*car_helper->input_dt*x1;
    x2(0) = u(1) * cos((double)xtmp(2));
    x2(1) = u(1) * sin((double)xtmp(2));
    x2(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
    xtmp = x + 0.5*car_helper->input_dt*x2;
    x3(0) = u(1) * cos((double)xtmp(2));
    x3(1) = u(1) * sin((double)xtmp(2));
    x3(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
    xtmp = x + car_helper->input_dt*x3;
    x4(0) = u(1) * cos((double)xtmp(2));
    x4(1) = u(1) * sin((double)xtmp(2));
    x4(2) = u(1) * tan((double)u(0)) / car_helper->carlen;

    new_x = x + car_helper->input_dt/6.0*(x1+2.0*(x2+x3)+x4);
    new_x.tail<4>() = x.tail<4>();
    new_x.head<3>() += 0.01 * m.head<3>();

    new_x.tail<4>() += 0.01 * m.tail<4>();

    double umag = u.squaredNorm();
    //return new_x + 0.1*umag*m;
    return new_x;// + 0.01*m;
  }

  CarObserveFunc::CarObserveFunc() : ObserveFunc<StateT, ObserveT, ObserveNoiseT>() {}

  CarObserveFunc::CarObserveFunc(BSPProblemHelperBasePtr helper) :
              ObserveFunc<StateT, ObserveT, ObserveNoiseT>(helper), car_helper(boost::static_pointer_cast<CarBSPProblemHelper>(helper)) {}

  ObserveT CarObserveFunc::operator()(const StateT& x, const ObserveNoiseT& n) const {
    ObserveT ret(observe_dim);
    Vector2d car_pos = x.head<2>(), l1 = x.middleRows<2>(3), l2 = x.middleRows<2>(5);
    double car_angle = x(2);

    ret(0) = (car_pos - l1).norm();
    ret(1) = atan2(car_pos.y() - l1.y(), car_pos.x() - l1.x()) - car_angle;
    ret(2) = (car_pos - l2).norm();
    ret(3) = atan2(car_pos.y() - l2.y(), car_pos.x() - l2.x()) - car_angle;

    //ret(0) = (car_pos(0) - l1(0));
    //ret(1) = (car_pos(1) - l1(1));
    //ret(2) = (car_pos(0) - l2(0));
    //ret(3) = (car_pos(1) - l2(1));

    return ret + 0.01 * n;
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

  ObserveT CarObserveFunc::observation_masks(const StateT& x, double approx_factor) const {
    ObserveT ret(observe_dim);
    Vector2d dists;
    sgndist(x, &dists);
    double tol = 0;
    if (approx_factor < 0) {
      ret(0) = ret(1) = dists(0) <= 0 ? 1 : 0;
      ret(2) = ret(3) = dists(1) <= 0 ? 1 : 0;
    } else {
      ret(0) = ret(1) = 1. - sigmoid(approx_factor * (dists(0) + tol));
      ret(2) = ret(3) = 1. - sigmoid(approx_factor * (dists(1) + tol));
    }

    return ret;
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
    QPen cvx_cov_pen(Qt::blue, 3, Qt::SolidLine);
    //QPen path_pen(Qt::blue, 3, Qt::SolidLine);
    QPen path_pen(Qt::red, 3, Qt::SolidLine);
    QPen pos_pen(Qt::blue, 4, Qt::SolidLine);
    QPen fov_pen(QColor(255, 255, 255, 10), 1, Qt::SolidLine);
    QBrush fov_brush(QColor(255, 255, 255, 30));
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::HighQualityAntialiasing);
    //painter.setPen(cvx_cov_pen);
    //for (int i = 0; i < states_opt.size(); ++i) {
    //  draw_ellipse(states_opt[i].head<2>(), sigmas_opt[i].topLeftCorner(2,2), painter, 0.5);
    //}
    painter.setPen(fov_pen);
    {
      QBrush prev_brush = painter.brush();
      painter.setBrush(fov_brush);
      Beam2D fov;
      for (int i = 0; i < (int)states_opt.size()-1; ++i) {
        fov = robot_fov(states_opt[i](0), states_opt[i](1), states_opt[i](2), car_helper->car_camera_depth, car_helper->car_camera_span_angle);
        //cout << "fov: " << fov << endl;
        draw_beam_2d(fov, painter);
      }
      painter.setBrush(prev_brush);
    }

    painter.setPen(path_pen);
    //for (int i = 0; i < states_opt.size() - 1; ++i) {
    //  draw_line(states_opt[i](0), states_opt[i](1), states_opt[i+1](0), states_opt[i+1](1), painter);
    //}

    for (int i = 0; i < (int)states_waypoints.size() - 1; ++i) {
      draw_line(states_waypoints[i](0), states_waypoints[i](1), states_waypoints[i+1](0), states_waypoints[i+1](1), painter);
    }

    //vector<Vector4d>& edges = car_helper->rrt_edges;
    //for (int i = 0; i < edges.size(); ++i) {
    //  draw_line(edges[i](0), edges[i](1), edges[i](2), edges[i](3), painter);
    //}
    //painter.setPen(pos_pen);
    //for (int i = 0; i < states_opt.size(); ++i) {
    //  draw_point(states_opt[i](0), states_opt[i](1), painter);
    ////for (int i = 0; i < samples.size(); ++i) {
    ////	draw_point(samples[i](0), samples[i](1), painter);
    //}

    QPen sensor_cov_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_cov_pen);
    //for (int i = 0; i < states_opt.size(); ++i) {
      draw_ellipse(states_actual.back().middleRows<2>(3), sigmas_actual.back().block<2, 2>(3, 3), painter, 0.5);
    //}
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < (int) states_actual.size() - 1; ++i) {
      draw_line(states_actual[i](3), states_actual[i](4), states_actual[i+1](3), states_actual[i+1](4), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_actual.size(); ++i) {
      draw_point(states_actual[i](3), states_actual[i](4), painter);
    }

    painter.setPen(sensor_cov_pen);
    //for (int i = 0; i < states_opt.size(); ++i) {
    //  draw_ellipse(states_opt[i].middleRows<2>(5), sigmas_opt[i].block<2, 2>(5, 5), painter);
    //}
      draw_ellipse(states_actual.back().middleRows<2>(5), sigmas_actual.back().block<2, 2>(5, 5), painter, 0.5);
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < (int) states_actual.size() - 1; ++i) {
      draw_line(states_actual[i](5), states_actual[i](6), states_actual[i+1](5), states_actual[i+1](6), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < states_actual.size(); ++i) {
      draw_point(states_actual[i](5), states_actual[i](6), painter);
    }

    // draw beliefs computed using belief dynamics
    QPen cvx_cov_pen2(QColor(255,215,0), 3, Qt::SolidLine);
    painter.setPen(cvx_cov_pen2);
    draw_ellipse(states_actual[0].head<2>(), sigmas_actual[0].topLeftCorner(2,2), painter, 0.5);
    draw_ellipse(states_actual.back().head<2>(), sigmas_actual.back().topLeftCorner(2,2), painter, 0.5);
    //for (int i = 0; i < states_actual.size(); ++i) {
    //  draw_ellipse(states_actual[i].head<2>(), sigmas_actual[i].topLeftCorner(2,2), painter);
    //}
    //QPen pos_pen2(Qt::red, 8, Qt::SolidLine);
    QPen pos_pen2(QColor(255,215,0), 8, Qt::SolidLine);
    painter.setPen(pos_pen2);
    for (int i = 0; i < states_actual.size(); ++i) {
      draw_point(states_actual[i](0), states_actual[i](1), painter);
    }

    // draw goal
    QPen goal_pen(QColor(192,192,192,100), 2, Qt::SolidLine);
    painter.setPen(goal_pen);
    QBrush goal_brush(QColor(192, 192, 192, 100));
    QBrush prev_brush = painter.brush();
    painter.setBrush(goal_brush);
    Vector2d g = car_helper->goal.head<2>();
    draw_ellipse(g, Matrix2d::Identity()*car_helper->goaleps, painter);
    painter.setBrush(prev_brush);
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
    states_waypoints.clear();

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

        StateT waypoint = cur_state;
        states_waypoints.push_back(waypoint);
        ControlT u = (ControlT) getVec(*xvec, car_helper->control_vars.row(i));
        double dt = 0.2*car_helper->input_dt;
        for(int j = 0; j < 5; ++j) {
          StateT xtmp(state_dim), x1(state_dim), x2(state_dim), x3(state_dim), x4(state_dim);
          xtmp = waypoint;
          x1(0) = u(1) * cos((double)xtmp(2));
          x1(1) = u(1) * sin((double)xtmp(2));
          x1(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + 0.5*dt*x1;
          x2(0) = u(1) * cos((double)xtmp(2));
          x2(1) = u(1) * sin((double)xtmp(2));
          x2(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + 0.5*dt*x2;
          x3(0) = u(1) * cos((double)xtmp(2));
          x3(1) = u(1) * sin((double)xtmp(2));
          x3(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + dt*x3;
          x4(0) = u(1) * cos((double)xtmp(2));
          x4(1) = u(1) * sin((double)xtmp(2));
          x4(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          waypoint = waypoint + dt/6.0*(x1+2.0*(x2+x3)+x4);
          states_waypoints.push_back(waypoint);
        }
      }
    }
    states_opt = new_states_opt; states_actual = new_states_actual;
    sigmas_opt = new_sigmas_opt; sigmas_actual = new_sigmas_actual;

    for(int i = 0; i < 1000; ++i) {
      samples.push_back(sample_gaussian(car_helper->start, car_helper->start_sigma, car_helper->noise_level));
    }

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
    for (int i = 0; i < (int)simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i+1](0), simulated_positions[i+1](1), painter);
    //for (int i = 0; i < (int)waypoints.size() - 1; ++i) {
      //draw_line(waypoints[i](0), waypoints[i](1), waypoints[i+1](0), waypoints[i+1](1), painter);
    }

    {
      QBrush prev_brush = painter.brush();
      painter.setPen(fov_pen);
      painter.setBrush(fov_brush);
      for (int i = 0; i < (int)simulated_positions.size(); ++i) {
        Beam2D fov = robot_fov(simulated_positions[i](0), simulated_positions[i](1), simulated_positions[i](2), car_helper->car_camera_depth, car_helper->car_camera_span_angle);
        draw_beam_2d(fov, painter);
      }
      painter.setBrush(prev_brush);
    }

    QPen sensor_path_pen(Qt::green, 3, Qt::SolidLine);
    QPen sensor_pos_pen(Qt::green, 8, Qt::SolidLine);
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < (int)simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](3), simulated_positions[i](4), simulated_positions[i+1](3), simulated_positions[i+1](4), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < (int)simulated_positions.size(); ++i) {
      draw_point(simulated_positions[i](3), simulated_positions[i](4), painter);
    }
    painter.setPen(sensor_path_pen);
    for (int i = 0; i < (int)simulated_positions.size() - 1; ++i) {
      draw_line(simulated_positions[i](5), simulated_positions[i](6), simulated_positions[i+1](5), simulated_positions[i+1](6), painter);
    }
    painter.setPen(sensor_pos_pen);
    for (int i = 0; i < (int)simulated_positions.size(); ++i) {
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
    waypoints.clear();

    for (int i = 0; i <= car_helper->T; ++i) {
      StateT cur_state;
      VarianceT cur_sigma;
      car_helper->belief_func->extract_state(cur_belief, &cur_state);
      car_helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      new_states.push_back(cur_state);
      new_sigmas.push_back(cur_sigma);
      if (i < car_helper->T) {
        cur_belief = car_helper->belief_func->call(cur_belief, (ControlT) getVec(*xvec, car_helper->control_vars.row(i)));
        /*
        StateT waypoint = cur_state;
        waypoints.push_back(waypoint);
        ControlT u = (ControlT) getVec(*xvec, car_helper->control_vars.row(i));
        double dt = 0.2*car_helper->input_dt;
        for(int j = 0; j < 5; ++j) {
          StateT xtmp(state_dim), x1(state_dim), x2(state_dim), x3(state_dim), x4(state_dim);
          xtmp = waypoint;
          x1(0) = u(1) * cos((double)xtmp(2));
          x1(1) = u(1) * sin((double)xtmp(2));
          x1(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + 0.5*dt*x1;
          x2(0) = u(1) * cos((double)xtmp(2));
          x2(1) = u(1) * sin((double)xtmp(2));
          x2(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + 0.5*dt*x2;
          x3(0) = u(1) * cos((double)xtmp(2));
          x3(1) = u(1) * sin((double)xtmp(2));
          x3(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          xtmp = waypoint + dt*x3;
          x4(0) = u(1) * cos((double)xtmp(2));
          x4(1) = u(1) * sin((double)xtmp(2));
          x4(2) = u(1) * tan((double)u(0)) / car_helper->carlen;
          waypoint = waypoint + dt/6.0*(x1+2.0*(x2+x3)+x4);
          waypoints.push_back(waypoint);
        }
        */
      }
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

    /* enum Method { StateSpace = 0, ContinuousBeliefSpace = 1, DiscontinuousBeliefSpace = 2}; */
    int method = 0;

    double noise_level = 1;

    double start_vec_array[] = {-5, 2, -PI*0.5, 0, 2, 0, -2};
    double goal_vec_array[] = {-5, -2, 0, 0, 0, 0, 0};

    vector<double> start_vec(start_vec_array, end(start_vec_array));
    vector<double> goal_vec(goal_vec_array, end(goal_vec_array));

    {
      Config config;
      config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
      config.add(new Parameter<int>("method", &method, "method"));
      config.add(new Parameter<double>("noise_level", &noise_level, "noise_level"));
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
    planner->method = method;
    planner->T = 30;
    planner->noise_level = noise_level;
    planner->initialize();
    

    boost::shared_ptr<CarSimulationPlotter> sim_plotter;
    boost::shared_ptr<CarOptPlotter> opt_plotter;

    if (plotting) {
      double x_min = -11.25, x_max = 6.25, y_min = -5, y_max = 5;
      sim_plotter.reset(create_plotter<CarSimulationPlotter>(x_min, x_max, y_min, y_max, planner->helper));
      sim_plotter->show();
      if (method != 0) {
        opt_plotter.reset(create_plotter<CarOptPlotter>(x_min, x_max, y_min, y_max, planner->helper));
        opt_plotter->show();
      }
    }

    boost::function<void(OptProb*, DblVec&)> opt_callback;
    if (plotting && method != 0) {
      opt_callback = boost::bind(&CarOptimizerTask::stage_plot_callback, this, opt_plotter, _1, _2);
    }

    //cout << "start solving" << endl;

    while (!planner->finished()) {
      planner->solve(opt_callback);
      cout << planner->helper->belief_func->approx_factor << endl;
      //planner->simulate_executions(planner->helper->T);
      planner->simulate_executions(1);
      if (plotting) {
        cout << "plotting now" << endl;
        emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions, true);
        //sim_plotter->update_plot_data(&planner->result, &planner->simulated_positions);
      }
    }

    cout << "final position: " << planner->simulated_positions.back().head<2>().transpose() << endl;
    cout << "goal position: " << planner->goal.head<2>().transpose() << endl;
    cout << "final estimated position: " << planner->start.head<2>().transpose() << endl;
    cout << "final covariance: " << endl << planner->start_sigma.topLeftCorner<2, 2>() << endl;
    if (plotting) {
      emit_plot_message(sim_plotter, &planner->result, &planner->simulated_positions);
      emit finished_signal();
    }
  }
}

using namespace CarBSP;

int main(int argc, char *argv[]) {
  seed_random();
  //srand(static_cast<unsigned int>(std::time(0)));
  bool plotting = true;
  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }
  QApplication app(argc, argv);
  CarOptimizerTask* task = new CarOptimizerTask(argc, argv, &app);
  if (plotting) {
    QTimer::singleShot(0, task, SLOT(run_slot()));
    QObject::connect(task, SIGNAL(finished_signal()), &app, SLOT(quit()));
    return app.exec();
  } else {
    task->run();
    return 0;
	}
}
