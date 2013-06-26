#include "bsp_optimizer_task.hpp"

namespace BSP {

  BSPOptimizerTask::BSPOptimizerTask(QObject *parent) : QObject(parent) {}

  BSPOptimizerTask::BSPOptimizerTask(int argc, char **argv, QObject *parent) : QObject(parent), argc(argc), argv(argv) {}

  void BSPOptimizerTask::run_slot() {
    run();
  }

  void BSPOptimizerTask::proceed_slot() {
    emit proceed_signal();
  }

  void BSPOptimizerTask::emit_plot_message(OptProb* prob, DblVec& xvec) {
    QEventLoop loop;
    connect(this, SIGNAL(proceed_signal()), &loop, SLOT(quit()));
    emit replot_signal(prob, xvec);
    loop.exec();
  }
}
