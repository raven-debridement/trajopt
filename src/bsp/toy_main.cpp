#include "bsp.hpp"
#include "toy.hpp"
#include "utils/logging.hpp"
#include "toy_optimizer_task.hpp"
#include "toy_plotter.hpp"
#include <QApplication>
#include <QtCore>
#include <qwt_plot_canvas.h>

#define CUSTOM_PREFIX "\x1b[32m[CUSTOM] "
#define LOG_CUSTOM(msg, ...) {printf(CUSTOM_PREFIX); printf(msg, ##__VA_ARGS__); printf(LOG_SUFFIX);}

using namespace BSP;
using namespace ToyBSP;

int main(int argc, char *argv[]) {

  QApplication app(argc, argv);

  bool plotting = true;
  {
    Config config;
    config.add(new Parameter<bool>("plotting", &plotting, "plotting"));
    CommandParser parser(config);
    parser.read(argc, argv, true);
  }

  ToyOptimizerTask* task = new ToyOptimizerTask(argc, argv, &app, plotting);
  ToyPlotter plotter;

  if (plotting) {
    QObject::connect(task, SIGNAL(plot(DblVec&, ToyBSPProblemHelperPtr)), &plotter, SLOT(plot(DblVec&, ToyBSPProblemHelperPtr)));
    QObject::connect(&plotter, SIGNAL(finished()), &app, SLOT(quit()));
    QObject::connect(&plotter, SIGNAL(proceed_signal()), task, SLOT(proceed_slot()));
    plotter.show();
  } else {
    QObject::connect(task, SIGNAL(finished()), &app, SLOT(quit()));
  }
  QTimer::singleShot(0, task, SLOT(run()));

  return app.exec();
}
