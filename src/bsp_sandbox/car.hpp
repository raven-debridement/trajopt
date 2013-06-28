#pragma once

#include "bsp/bsp.hpp"
#include <QtCore>
#include <QImage>
#include <qevent.h>

using namespace BSP;

namespace CarBSP {

// This will generate a bunch of types like StateT, ControlT, etc.
BSP_TYPEDEFS(
		4, // state_dim
		4, // state_noise_dim
		2, // control_dim
		2, // observe_dim
		2, // observe_noise_dim
		10, // sigma_dof
		14 // belief_dim
);

// state: { x, y, angle, velocity }
// control: { theta, acceleration }

class CarBSPProblemHelper;
typedef boost::shared_ptr<CarBSPProblemHelper> CarBSPProblemHelperPtr;

class CarStateFunc : public StateFunc<StateT, ControlT, StateNoiseT> {
public:
	typedef boost::shared_ptr<CarStateFunc> Ptr;
	CarBSPProblemHelperPtr car_helper;

	CarStateFunc();
	CarStateFunc(BSPProblemHelperBasePtr helper);
	virtual StateT operator()(const StateT& x, const ControlT& u, const StateNoiseT& m) const ;
};

class CarObserveFunc : public ObserveFunc<StateT, ObserveT, ObserveNoiseT> {
public:
	typedef boost::shared_ptr<CarObserveFunc> Ptr;
	CarBSPProblemHelperPtr car_helper;

	CarObserveFunc();
	CarObserveFunc(BSPProblemHelperBasePtr helper);
	virtual ObserveT operator()(const StateT& x, const ObserveNoiseT& n) const;
};

class CarBeliefFunc : public BeliefFunc<CarStateFunc, CarObserveFunc, BeliefT> {
public:
	typedef boost::shared_ptr<CarBeliefFunc> Ptr;
	CarBSPProblemHelperPtr car_helper;
	double alpha;
	double tol;

	CarBeliefFunc();
	CarBeliefFunc(BSPProblemHelperBasePtr helper, StateFuncPtr f, ObserveFuncPtr h);
	bool sgndist(const Vector2d& x, Vector2d* dists) const;
	virtual ObserveStateGradT sensor_constrained_observe_state_gradient(const ObserveStateGradT& H, const StateT& x) const;
	virtual ObserveNoiseGradT sensor_constrained_observe_noise_gradient(const ObserveNoiseGradT& N, const StateT& x) const;
	virtual VarianceT sensor_constrained_variance_reduction(const VarianceT& reduction, const StateT& x) const;
	virtual ObserveMatT compute_gamma(const StateT& x) const;
	virtual ObserveMatT compute_inverse_gamma(const StateT& x) const;
};

class CarBSPProblemHelper : public BSPProblemHelper<CarBeliefFunc> {
public:
	struct RRTNode {
		StateT x;
		ControlT u;
		int bp;
	};

	typedef typename BeliefConstraint<CarBeliefFunc>::Ptr BeliefConstraintPtr;
	double input_dt;
	double carlen;
	vector<ControlT> init_controls;

	virtual void RRTplan();
	virtual void init_control_values(vector<ControlT>* output_init_controls) const;
	virtual void add_goal_constraint(OptProb& prob);
	CarBSPProblemHelper();
};

class CarPlotter : public BSPQtPlotter, public ProblemState {
	Q_OBJECT
public:
	CarPlotter(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper, QWidget* parent=NULL);
public slots:
virtual void update_plot_data(OptProb*, DblVec& x);
protected:
double old_alpha, cur_alpha;
CarBSPProblemHelperPtr car_helper;
QImage distmap;
vector<VectorXd> states;
vector<MatrixXd> sigmas;
virtual void paintEvent(QPaintEvent*);
};

class CarOptimizerTask : public BSPOptimizerTask {
	Q_OBJECT
public:
	CarOptimizerTask(QObject* parent=NULL);
	CarOptimizerTask(int argc, char **argv, QObject* parent=NULL);
	virtual void emit_plot_message(OptProb* prob, DblVec& xvec);
	virtual void run();
};
}

template CarBSP::CarPlotter* BSPOptimizerTask::create_plotter<CarBSP::CarPlotter>(double x_min, double x_max, double y_min, double y_max, BSPProblemHelperBasePtr helper);
