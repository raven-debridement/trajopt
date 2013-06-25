#include "toy_plotter.hpp"


namespace ToyBSP {

  extern QWaitCondition optimizer_proceed;

  class CoordinateScaler {
  public:
    CoordinateScaler(double x_min, double x_max, double y_min, double y_max, double width, double height) :
      x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), width(width), height(height) {}

    double scale_x_length(double x) {
      return x / (x_max - x_min) * width;
    }

    double scale_y_length(double y) {
      return y / (y_max - y_min) * height;
    }

    double scale_x(double x) {
      return scale_x_length(x - x_min);
    }

    double scale_y(double y) {
      return height - scale_y_length(y - y_min);
    }

    double unscale_x(double x) {
      return (x / width) * (x_max - x_min) + x_min;
    }

    double unscale_y(double y) {
      return ((height - y) / height) * (y_max - y_min) + y_min;
    }

  protected:
    double x_min, x_max, y_min, y_max, width, height;
  };

  ToyPlotter::ToyPlotter(QWidget* parent) : QWidget(parent), old_alpha(-1), cur_alpha(-1), distmap(400, 400, QImage::Format_RGB32) {
    QPalette pal = palette();
    pal.setColor(backgroundRole(), Qt::white);
    setPalette(pal);
  }

  void ToyPlotter::paintEvent(QPaintEvent* ) {
    if (!helper) {
      cout << "helper undefined\n";
      return;
    }
    double x_min = -7, x_max = 2, y_min = -1, y_max = 3;
    double ratio = height() * 1.0 / width();
    CoordinateScaler scaler(x_min, x_max, y_min, y_max, width(), height());
    QPainter painter(this);
    
    if (cur_alpha != old_alpha || distmap.height() != height() || distmap.width() != width()) { // replot distmap
      distmap = QImage(width(), height(), QImage::Format_RGB32);
      for (int j = 0; j < height(); ++j) {
        QRgb *line = (QRgb*) distmap.scanLine(j);
        for (int i = 0; i < width(); ++i) {
          double x = scaler.unscale_x(i),
                 y = scaler.unscale_y(j);
          StateT dists(helper->state_dim);
          StateT state(helper->state_dim); state << x, y;
          helper->belief_func->sgndist(state, &dists);
          double grayscale = fmax(1./(1. + exp(helper->belief_func->alpha*dists(0))),
                                  1./(1. + exp(helper->belief_func->alpha*dists(1))));
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
    for (int i = 0; i < ellipse_params.size(); ++i) {
      double cx = scaler.scale_x(ellipse_params[i](0)),
             cy = scaler.scale_y(ellipse_params[i](1)),
             mx = scaler.scale_x_length(ellipse_params[i](2)),
             my = scaler.scale_y_length(ellipse_params[i](3)),
             sxy = ellipse_params[i](4),
             sxy_diff = ellipse_params[i](5);
      double theta = 0.5 * atan2(2*sxy, sxy_diff * ratio);
      painter.save();
      painter.translate(cx-mx, cy-my);
      painter.rotate(theta);
      painter.drawEllipse(0, 0, 2*mx, 2*my);  
      painter.restore();
    }

    painter.setPen(path_pen);

    for (int i = 0; i < ellipse_params.size() - 1; ++i) {
      double curx = scaler.scale_x(ellipse_params[i](0)),
             cury = scaler.scale_y(ellipse_params[i](1)),
             nextx = scaler.scale_x(ellipse_params[i+1](0)),
             nexty = scaler.scale_y(ellipse_params[i+1](1));
      painter.drawLine(curx, cury, nextx, nexty);
    }

    painter.setPen(pos_pen);

    for (int i = 0; i < ellipse_params.size(); ++i) {
      double curx = scaler.scale_x(ellipse_params[i](0)),
             cury = scaler.scale_y(ellipse_params[i](1));
      painter.drawPoint(curx, cury);
    }

    
  }

  void ToyPlotter::plot(DblVec& xvec, ToyBSPProblemHelperPtr helper) {
    vector<VectorXd> new_ellipse_params;
    old_alpha = cur_alpha;
    cur_alpha = helper->belief_func->alpha;
    this->helper = helper;
    
    BeliefT cur_belief;
    
    helper->belief_func->compose_belief(helper->start, helper->start_sigma, &cur_belief);

    for (int i = 0; i <= helper->T; ++i) {
      VarianceT cur_sigma;
      helper->belief_func->extract_sigma(cur_belief, &cur_sigma);
      EigenSolver<VarianceT> es(cur_sigma);
      StateT eigenvals = es.eigenvalues().real();
      double scale_factor = 0.5;//2.4477;
      double sx = cur_sigma(0, 0), sy = cur_sigma(1, 1);
      double max_val = eigenvals(0), min_val = eigenvals(1);
      if (min_val > max_val) {
        double tmp = max_val; max_val = min_val; min_val = tmp;
      }
      double sxy = cur_sigma(0, 1);
      double sxy_diff = sx - sy;
      //double theta = 0.5 * atan2(2*cur_sigma(0,1), sx - sy) * PI;
      double mx, my;
      if (sx >= sy) {
        mx = sqrt(max_val)*scale_factor;
        my = sqrt(min_val)*scale_factor;
      } else {
        mx = sqrt(min_val)*scale_factor;
        my = sqrt(max_val)*scale_factor;
      }
      double cx = cur_belief(0), cy = cur_belief(1);
      if (i < helper->T) cur_belief = helper->belief_func->call(cur_belief, (ControlT) getVec(xvec, helper->control_vars.row(i)));
      VectorXd ellipse_param(6); ellipse_param << cx, cy, mx, my, sxy, sxy_diff;
      new_ellipse_params.push_back(ellipse_param);
    }
    ellipse_params = new_ellipse_params;
    this->repaint();
  }

  void ToyPlotter::keyPressEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_P) {
      emit(proceed_signal());
    }
  }

  void ToyPlotter::closeEvent(QCloseEvent*) {
    emit finished();
  }

  void ToyPlotter::showEvent(QShowEvent*) {
    emit inited();
  }
}
