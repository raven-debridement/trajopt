#pragma once

#include <QtCore>
#include "common.hpp"
#include "bsp_problem_helper_base.hpp"

namespace BSP {
  class BSPOptimizerTask : public QObject {
    Q_OBJECT
  public:
    BSPOptimizerTask(QObject* parent=NULL);
    BSPOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void emit_plot_message(OptProb* prob, DblVec& xvec);
    virtual void run() = 0;
    int argc;
    char **argv;
  public slots:
    void run_slot();
    void proceed_slot();
  signals:
    void finished_signal();
    void proceed_signal();
    void replot_signal(OptProb* prob, DblVec& xvec);
  protected:
    template<class PlotterT> PlotterT* create_plotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper) {
      PlotterT* plotter = new PlotterT(x_min, x_max, y_min, y_max, helper, NULL);
      QObject::connect(this, SIGNAL(replot_signal(OptProb*, DblVec&)), plotter, SLOT(update_plot_data(OptProb*, DblVec&)));
      QObject::connect(plotter, SIGNAL(finished_signal()), parent(), SLOT(quit()));
      QObject::connect(plotter, SIGNAL(proceed_signal()), this, SLOT(proceed_slot()));
      return plotter;
    }
  };
}
