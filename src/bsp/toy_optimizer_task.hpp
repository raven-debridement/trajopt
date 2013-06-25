#pragma once

#include "bsp.hpp"
#include "toy.hpp"
#include "utils/logging.hpp"
#define CUSTOM_PREFIX "\x1b[32m[CUSTOM] "
#define LOG_CUSTOM(msg, ...) {printf(CUSTOM_PREFIX); printf(msg, ##__VA_ARGS__); printf(LOG_SUFFIX);}

#include <QtCore>
#include <qevent.h>
using namespace BSP;

namespace ToyBSP {

  class ToyOptimizerTask : public QObject {
    Q_OBJECT
  public:
    ToyOptimizerTask(QObject *parent=0);
    ToyOptimizerTask(int argc, char **argv, QObject *parent=0, bool plotting=false);
    void emit_plot_message(OptProb*, DblVec& xvec, ToyBSPProblemHelperPtr helper);
  public slots:
    void run();
    void proceed_slot();
  signals:
    void finished();
    void proceed_signal();
    void plot(DblVec& xvec, ToyBSPProblemHelperPtr helper);
  protected:
    int argc;
    char **argv;
    bool plotting;
  };
}
