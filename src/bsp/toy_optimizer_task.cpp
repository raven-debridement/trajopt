#include "toy_optimizer_task.hpp"

namespace ToyBSP {

  QMutex plot_mutex;
  QWaitCondition optimizer_proceed;
  ToyOptimizerTask::ToyOptimizerTask(QObject *parent) : QObject(parent), plotting(false) {}

  ToyOptimizerTask::ToyOptimizerTask(int argc, char **argv, QObject *parent, bool plotting) : QObject(parent), argc(argc), argv(argv), plotting(plotting) {}

  void ToyOptimizerTask::proceed_slot() {
    emit proceed_signal();
  }

  void ToyOptimizerTask::emit_plot_message(OptProb*, DblVec& xvec, ToyBSPProblemHelperPtr helper) {
    QEventLoop loop;
    connect(this, SIGNAL(proceed_signal()), &loop, SLOT(quit()));
    emit plot(xvec, helper);
    loop.exec();
  }

  void ToyOptimizerTask::run() {
    int T = 20;
    int n_dof = 2;

    double start_vec_array[] = {-5, 2};
    double goal_vec_array[] = {-5, 0};

    vector<double> start_vec(start_vec_array, end(start_vec_array));
    vector<double> goal_vec(goal_vec_array, end(goal_vec_array));

    {
      Config config;
      config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
      config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
      CommandParser parser(config);
      parser.read(argc, argv, true);
    }

    Vector2d start = toVectorXd(start_vec);
    Vector2d goal = toVectorXd(goal_vec);
    Matrix2d start_sigma = Matrix2d::Identity();

    OptProbPtr prob(new OptProb());

    ToyBSPProblemHelperPtr helper(new ToyBSPProblemHelper());
    helper->start = start;
    helper->goal = goal;
    helper->start_sigma = start_sigma;
    helper->T = T;
    helper->configure_problem(*prob);

    BSPTrustRegionSQP opt(prob);
    opt.max_iter_ = 500;
    opt.merit_error_coeff_ = 500;
    opt.trust_shrink_ratio_ = 0.5;
    opt.trust_expand_ratio_ = 1.25;
    opt.min_trust_box_size_ = 1e-3;
    opt.min_approx_improve_ = 1e-2;
    opt.min_approx_improve_frac_ = 1e-4;
    opt.improve_ratio_threshold_ = 0.2;
    opt.trust_box_size_ = 1;

    helper->configure_optimizer(*prob, opt);

    if (plotting) {
      opt.addCallback(boost::bind(&ToyOptimizerTask::emit_plot_message, boost::ref(this), _1, _2, helper));
    }

    opt.optimize();

    emit finished();
  }
}
