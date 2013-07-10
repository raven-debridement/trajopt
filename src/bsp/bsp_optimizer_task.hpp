#pragma once

#include <QtCore>
#include "common.hpp"
#include "bsp_problem_helper_base.hpp"

namespace BSP {
  class BSPOptimizerTask : public QObject {
    Q_OBJECT
  public:
    int argc;
    char **argv;
    BSPOptimizerTask(QObject* parent=NULL);
    BSPOptimizerTask(int argc, char **argv, QObject* parent=NULL);
    virtual void run() = 0;

    void wait_to_proceed(boost::function<void()> f, bool wait=true) {
      if (wait) {
        QEventLoop loop;
        connect(this, SIGNAL(proceed_signal()), &loop, SLOT(quit()));
        f();
        loop.exec();
      } else {
        f();
      }
    }

    template<class PlotterT>
    void emit_plot_message(boost::shared_ptr<PlotterT> plotter, void* d1, bool wait=true) {
      wait_to_proceed(boost::bind(&PlotterT::update_plot_data, plotter, d1), wait);
    }

    template<class PlotterT>
    void emit_plot_message(boost::shared_ptr<PlotterT> plotter, void* d1, void* d2, bool wait=true) {
      wait_to_proceed(boost::bind(&PlotterT::update_plot_data, plotter, d1, d2), wait);
    }

    template<class PlotterT>
    void emit_plot_message(boost::shared_ptr<PlotterT> plotter, void* d1, void* d2, void* d3, bool wait=true) {
      wait_to_proceed(boost::bind(&PlotterT::update_plot_data, plotter, d1, d2, d3), wait);
    }

  public slots:
    void run_slot();
    void proceed_slot();
  signals:
    void finished_signal();
    void proceed_signal();
  protected:
    template<class PlotterT> PlotterT* create_plotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper) {
      PlotterT* plotter = new PlotterT(x_min, x_max, y_min, y_max, helper, NULL);
      QObject::connect(plotter, SIGNAL(finished_signal()), parent(), SLOT(quit()));
      QObject::connect(plotter, SIGNAL(proceed_signal()), this, SLOT(proceed_slot()));
      return plotter;
    }
  };
}
