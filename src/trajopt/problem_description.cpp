#include "trajopt/problem_description.hpp"
#include "trajopt/common.hpp"
#include <boost/foreach.hpp>
#include "utils/logging.hpp"
#include "sco/expr_ops.hpp"
#include "trajopt/kinematic_constraints.hpp"
#include "trajopt/belief_constraints.hpp"
#include "trajopt/collision_avoidance.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/plot_callback.hpp"
#include "trajopt/rave_utils.hpp"
#include "utils/eigen_conversions.hpp"
#include "utils/eigen_slicing.hpp"
#include <boost/algorithm/string.hpp>
using namespace Json;
using namespace std;
using namespace OpenRAVE;
using namespace trajopt;
using namespace util;

namespace {


bool gRegisteredMakers = false;
void RegisterMakers() {

	CostInfo::RegisterMaker("pose", &PoseCostInfo::create);
	CostInfo::RegisterMaker("joint_pos", &JointPosCostInfo::create);
	CostInfo::RegisterMaker("joint_vel", &JointVelCostInfo::create);
	CostInfo::RegisterMaker("collision", &CollisionCostInfo::create);
	CostInfo::RegisterMaker("continuous_collision", &ContinuousCollisionCostInfo::create);

	CntInfo::RegisterMaker("joint", &JointConstraintInfo::create);
	CntInfo::RegisterMaker("pose", &PoseCntInfo::create);
	CntInfo::RegisterMaker("collision", &CollisionCntInfo::create);
	CntInfo::RegisterMaker("continuous_collision", &ContinuousCollisionCntInfo::create);
	CntInfo::RegisterMaker("cart_vel", &CartVelCntInfo::create);

	// belief costs and controls
	CostInfo::RegisterMaker("control", &ControlCostInfo::create);
	CostInfo::RegisterMaker("covariance", &CovarianceCostInfo::create);

	CntInfo::RegisterMaker("control", &ControlCntInfo::create);

	gRegisteredMakers = true;
}

BeliefRobotAndDOFPtr RADFromName(const string& name, RobotBasePtr robot) {
	if (name == "active") {
		return BeliefRobotAndDOFPtr(new BeliefRobotAndDOF(robot, robot->GetActiveDOFIndices(), robot->GetAffineDOF(), robot->GetAffineRotationAxis()));
	}
	vector<int> dof_inds;
	int affinedofs = 0;
	Vector rotationaxis(0,0,1);
	vector<string> components;
	boost::split(components, name, boost::is_any_of("+"));
	for (int i=0; i < components.size(); ++i) {
		std::string& component = components[i];
		if (RobotBase::ManipulatorPtr manip = GetManipulatorByName(*robot, component)) {
			vector<int> inds = manip->GetArmIndices();
			dof_inds.insert(dof_inds.end(), inds.begin(), inds.end());
		}
		else if (component == "base") {
			affinedofs |= DOF_X | DOF_Y | DOF_RotationAxis;
		}
		else if (component == "base_point") {
			affinedofs |= DOF_X | DOF_Y;
		}
		else if (KinBody::JointPtr joint = robot->GetJoint(component)) {
			dof_inds.push_back(joint->GetDOFIndex());
		}
		else PRINT_AND_THROW( boost::format("error in reading manip description: %s must be a manipulator, link, or 'base'")%component );
	}
	return BeliefRobotAndDOFPtr(new BeliefRobotAndDOF(robot, dof_inds, affinedofs, rotationaxis));
}

BoolVec toMask(const VectorXd& x) {
	BoolVec out(x.size());
	for (int i=0; i < x.size(); ++i) out[i] = (x[i] > 0);
	return out;
}

bool allClose(const VectorXd& a, const VectorXd& b) {
	return (a-b).array().abs().maxCoeff() < 1e-4;
}

}

namespace Json { //funny thing with two-phase lookup

void fromJson(const Json::Value& v, Vector3d& x) {
	vector<double> vx;
	fromJsonArray(v, vx, 3);
	x = Vector3d(vx[0], vx[1], vx[2]);
}
void fromJson(const Json::Value& v, Vector4d& x) {
	vector<double> vx;
	fromJsonArray(v, vx, 4);
	x = Vector4d(vx[0], vx[1], vx[2], vx[3]);
}
template <class T>
inline void fromJson(const Json::Value& v, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>& m) {
	int nRows = v.size();
	if (nRows != 0) {
		int nCols = v[0].size();
		m.resize(nRows, nCols);
		for (int i=0; i<nRows; i++) {
			if (v[i].size() != nCols)
				PRINT_AND_THROW(boost::format("matrix with variable number of cols not supported. expected %i cols at the %ith row, but got %i cols\n")%nCols%i%v[i].size());
			std::vector<T> row;
			fromJsonArray(v[i], row, nCols);
			for (int j=0; j<nCols; j++) m(i,j) = row[j];
		}
	} else {
		m.resize(0,0);
	}
}

}

namespace trajopt {

TRAJOPT_API ProblemConstructionInfo* gPCI;

void BasicInfo::fromJson(const Json::Value& v) {
	childFromJson(v, start_fixed, "start_fixed", true);
	childFromJson(v, n_steps, "n_steps");
	childFromJson(v, manip, "manip");
	childFromJson(v, robot, "robot", string(""));
	childFromJson(v, dofs_fixed, "dofs_fixed", IntVec());
	childFromJson(v, belief_space, "belief_space", false);
}


////
void fromJson(const Json::Value& v, CostInfoPtr& cost) {
	string type;
	childFromJson(v, type, "type");
	cost = CostInfo::fromName(type);
	if (!cost) PRINT_AND_THROW( boost::format("failed to construct cost named %s")%type );
	cost->fromJson(v);
	childFromJson(v, cost->name, "name", type);
}
CostInfoPtr CostInfo::fromName(const string& type) {
	if (!gRegisteredMakers) RegisterMakers();
	if (name2maker.find(type) != name2maker.end()) {
		return (*name2maker[type])();
	}
	else {
		RAVELOG_ERROR("There is no cost of type%s\n", type.c_str());
		return CostInfoPtr();
	}
}
map<string, CostInfo::MakerFunc> CostInfo::name2maker;

void CostInfo::RegisterMaker(const std::string& type, MakerFunc f) {
	name2maker[type] = f;
}


////
//// almost copied


void fromJson(const Json::Value& v, CntInfoPtr& cnt) {
	string type;
	childFromJson(v, type, "type");
	LOG_DEBUG("reading constraint: %s", type.c_str());
	cnt = CntInfo::fromName(type);
	if (!cnt) PRINT_AND_THROW( boost::format("failed to construct constraint named %s")%type );
	cnt->fromJson(v);
	childFromJson(v, cnt->name, "name", type);
}
CntInfoPtr CntInfo::fromName(const string& type) {
	if (!gRegisteredMakers) RegisterMakers();
	if (name2maker.find(type) != name2maker.end()) {
		return (*name2maker[type])();
	}
	else {
		RAVELOG_ERROR("There is no constraint of type%s\n", type.c_str());
		return CntInfoPtr();
	}
}
map<string,CntInfo::MakerFunc> CntInfo::name2maker;

void InitInfo::fromJson(const Json::Value& v) {
	string type_str;
	childFromJson(v, type_str, "type");
	int n_steps = gPCI->basic_info.n_steps;
	int n_dof = gPCI->rad->GetDOF();
	bool belief_space = gPCI->basic_info.belief_space;
	MatrixXd rt_Sigma0;
	if (belief_space && (type_str == "stationary" || type_str == "given_traj")) {
		FAIL_IF_FALSE(v.isMember("initial_rt_sigma"));
		const Value& rt_Sigma0_v = v["initial_rt_sigma"];
		Json::fromJson(rt_Sigma0_v, rt_Sigma0);
		if (rt_Sigma0.rows()!=n_dof) PRINT_AND_THROW("initial square root of sigma has wrong number of rows");
		if (rt_Sigma0.cols()!=n_dof) PRINT_AND_THROW("initial square root of sigma has wrong number of cols");
	}

	if (type_str == "stationary") {
		if (!belief_space) {
			data = toVectorXd(gPCI->rad->GetDOFValues()).transpose().replicate(n_steps, 1);
		} else {
			data.resize(n_steps, B_DIM+U_DIM);
			VecBd theta;
			//gPCI->rad->composeBelief(toVectorXd(gPCI->rad->GetDOFValues()), rt_Sigma0, theta);
			VecXd x;
			const DblVec& xvec = gPCI->rad->GetDOFValues();
			for(int i = 0; i < xvec.size(); ++i) { x[i] = xvec[i]; }
			cout << x << endl;
			cout << rt_Sigma0 << endl;
			cout << theta << endl;

			gPCI->rad->composeBelief(x, rt_Sigma0, theta);
			for (int i=0; i < n_steps; i++) {
				data.block(i,0,1,B_DIM) = theta.transpose();
				if (i != (n_steps-1)) {
					VecUd u = VecUd::Zero();
					data.block(i,B_DIM,1,U_DIM) = u.transpose();
					theta = gPCI->rad->BeliefDynamics(theta,u);
				} else {
					data.block(i,B_DIM,1,U_DIM) = VecUd::Zero().transpose();
				}
			}
		}
	}
	else if (type_str == "given_traj") {
		FAIL_IF_FALSE(v.isMember("data"));
		const Value& data_v = v["data"];
		TrajArray x_data;
		Json::fromJson(data_v, x_data);
		if (x_data.rows()!=n_steps) PRINT_AND_THROW("x data has wrong number of rows");

		if (!belief_space) {
			if (x_data.cols()!=n_dof) PRINT_AND_THROW("x data has wrong number of cols");
			data = x_data;
		} else {
			if (x_data.cols() == B_DIM+U_DIM) {
				data = x_data;
			} else if (x_data.cols() == n_dof) {
				data.resize(n_steps, B_DIM+U_DIM);
				VecBd theta;
				gPCI->rad->composeBelief(x_data.row(0).transpose(), rt_Sigma0, theta);
				for (int i=0; i < n_steps; i++) {
					data.block(i,0,1,B_DIM) = theta.transpose();
					if (i != (n_steps-1)) {
						VecUd u = x_data.row(i+1).transpose() - x_data.row(i).transpose();
						data.block(i,B_DIM,1,U_DIM) = u.transpose();
						theta = gPCI->rad->BeliefDynamics(theta,u);
					} else {
						data.block(i,B_DIM,1,U_DIM) = VecUd::Zero().transpose();
					}
				}
			} else {
				PRINT_AND_THROW("x data has wrong number of cols");
			}
		}
	}
	else if (type_str == "straight_line") {
		FAIL_IF_FALSE(!belief_space);
		FAIL_IF_FALSE(v.isMember("endpoint"));
		DblVec endpoint;
		childFromJson(v, endpoint, "endpoint");
		if (endpoint.size() != n_dof) {
			PRINT_AND_THROW(boost::format("wrong number of dof values in initialization. expected %i got %j")%n_dof%endpoint.size());
		}
		data = TrajArray(n_steps, n_dof);
		DblVec start = gPCI->rad->GetDOFValues();
		for (int idof = 0; idof < n_dof; ++idof) {
			data.col(idof) = VectorXd::LinSpaced(n_steps, start[idof], endpoint[idof]);
		}
	}
}

void ProblemConstructionInfo::fromJson(const Value& v) {
	childFromJson(v, basic_info, "basic_info");

	RobotBasePtr robot = (basic_info.robot=="") ? GetRobot(*env) : GetRobotByName(*env, basic_info.robot);
	if (!robot) {
		PRINT_AND_THROW("couldn't get robot");
	}
	rad = RADFromName(basic_info.manip, robot);
	if (!rad) {
		PRINT_AND_THROW( boost::format("couldn't get manip %s")%basic_info.manip );
	}

	gPCI = this;
	if (v.isMember("costs")) fromJsonArray(v["costs"], cost_infos);
	if (v.isMember("constraints")) fromJsonArray(v["constraints"], cnt_infos);

	childFromJson(v, init_info, "init_info");
	gPCI = NULL;

}
void CntInfo::RegisterMaker(const std::string& type, MakerFunc f) {
	name2maker[type] = f;
}

TrajOptResult::TrajOptResult(OptResults& opt, TrajOptProb& prob) :
		  cost_vals(opt.cost_vals),
		  cnt_viols(opt.cnt_viols) {
	BOOST_FOREACH(const CostPtr& cost, prob.getCosts()) {
		cost_names.push_back(cost->name());
	}
	BOOST_FOREACH(const ConstraintPtr& cnt, prob.getConstraints()) {
		cnt_names.push_back(cnt->name());
	}
	traj = getTraj(opt.x, prob.GetVars());
}

Vector3d endEffectorPosition(BeliefRobotAndDOFPtr brad, VecXd dofs) {
	Vector3d eetrans;
	brad->ForwardKinematics(dofs, eetrans);

	double c6 = cos((double)dofs[6]);
	double s6 = sin((double)dofs[6]);

	eetrans[0] = c6*eetrans[0] - s6*eetrans[1];
	eetrans[0] = s6*eetrans[0] + c6*eetrans[1];
	eetrans[2] = eetrans[2] + 0.16;

	return eetrans;
}

void renderSigmaPts(BeliefRobotAndDOFPtr rad, const MatrixXd& sigma_pts, const osg::Vec4f& colorvec, vector<GraphHandlePtr>& handles, OSGViewerPtr viewer) {
	vector<KinBody::LinkPtr> links;
	vector<int> joint_inds;
	rad->GetAffectedLinks(links, true, joint_inds);

	// render sigma points
//	for (int j=0; j<sigma_pts.cols(); j++) {
//		rad->SetDOFValues(toDblVec(sigma_pts.col(j)));
//		handles.push_back(viewer->PlotKinBody(rad->GetRobot()));
////		if (j==0) SetColor(handles.back(), osg::Vec4f(0,0,1,1));
////		else //SetColor(handles.back(), osg::Vec4f(0,0,1,1));
//		SetColor(handles.back(), colorvec);
//	}

	// render convex hulls of sigma points
	vector<DblVec> dofvals(sigma_pts.cols());
	for (int i=0; i<sigma_pts.cols(); i++)
		dofvals[i] = toDblVec(sigma_pts.col(i));
	//cc->SetContactDistance(100);
	boost::shared_ptr<CollisionChecker> cc = CollisionChecker::GetOrCreate(*rad->GetRobot()->GetEnv());
	cc->PlotCastHull(*rad, links, dofvals, handles);
}

void renderTrajectory(BeliefRobotAndDOFPtr brad, OSGViewerPtr viewer, const TrajArray& traj, vector<GraphHandlePtr>& handles, int sample_rate = 1) {
//	viewer->GetEnv()->Load("/home/alex/rll/trajopt/data/barrett_sensor_wall.env.xml");
//	KinBodyPtr wall = viewer->GetEnv()->GetKinBody("wall");
//	handles.push_back(viewer->PlotKinBody(wall));
//	SetColor(handles.back(), osg::Vec4f(0.8,0.8,0,1));
//	SetTransparency(handles.back(), 0.1);
//	viewer->GetEnv()->Remove(wall);

	int n_skipped = sample_rate;
	Vector3d last_ee_pos = endEffectorPosition(brad, traj.block(0,0,1,X_DIM).transpose());
	for (int i=0; i<traj.rows(); i++) {
		if ((n_skipped < sample_rate) && ((endEffectorPosition(brad, traj.block(i,0,1,X_DIM).transpose())-last_ee_pos).norm()<0.1) && (i != (traj.rows()-1))) {
			n_skipped++;
			continue;
		}
		n_skipped = 0;
		last_ee_pos = endEffectorPosition(brad, traj.block(i,0,1,X_DIM).transpose());

		brad->SetDOFValues(toDblVec(traj.block(i,0,1,X_DIM).transpose()));
		handles.push_back(viewer->PlotKinBody(brad->GetRobot()));
		if (i==0 || i==(traj.rows()-1)) renderSigmaPts(brad, brad->sigmaPoints(traj.block(i,0,1,B_DIM).transpose()), osg::Vec4f(0,1,0,0.2), handles, viewer);
		viewer->Idle();
	}

	viewer->Idle();
	viewer->Idle();
	viewer->Idle();
	viewer->Idle();
}

bool isTrajectoryInCollision(CollisionCheckerPtr cc, TrajArray traj, BeliefRobotAndDOFPtr rad) {
	vector<Collision> collisions;
	cc->DiscreteCheckTrajectory(traj, rad, collisions);
	for (int i=0; i<collisions.size(); i++) {
		if (collisions[i].distance < 0) return true;
	}
	return false;
}

double calcTrajectoryCost(BeliefRobotAndDOFPtr brad, const TrajArray& traj) {
	VecBd theta;
	VecUd u;
	VecXd x;
	MatXXd rt_Sigma;
	MatXXd Sigma;
	MatXXd Q = MatXXd::Identity()*10;
	MatUUd R = MatUUd::Identity()*0.1;
	double total_cost = 0;
	for (int i=0; i<traj.rows(); i++) {
		theta = traj.block(i,0,1,B_DIM).transpose();
		u = traj.block(i,B_DIM,1,U_DIM).transpose();
		brad->decomposeBelief(theta, x, rt_Sigma);
		total_cost += ((MatXXd) (Q*rt_Sigma*rt_Sigma.transpose())).trace() + u.transpose()*R*u;
	}
	return total_cost;
}

VectorXd SimulateAndReplan(const Json::Value& root, OpenRAVE::EnvironmentBasePtr env, bool sigma_pts_scale, bool interactive) {
	ProblemConstructionInfo pci(env);
	pci.fromJson(root);

	TrajArray& traj = pci.init_info.data;
	TrajArray init_traj = traj;
	BeliefRobotAndDOFPtr brad = pci.rad;
	brad->SetSigmaPointsScale(sigma_pts_scale);
	int n_steps = pci.basic_info.n_steps;
	TrajArray plan_traj = TrajArray::Zero(n_steps, B_DIM+U_DIM);
	TrajArray exec_mpc_traj = TrajArray::Zero(n_steps, B_DIM+U_DIM);
	TrajArray exec_mpc_gt_traj = TrajArray::Zero(n_steps, X_DIM);
	TrajArray exec_open_traj = TrajArray::Zero(n_steps, B_DIM+U_DIM);
	TrajArray exec_open_gt_traj = TrajArray::Zero(n_steps, X_DIM);
//	TrajArray exec_mpc_traj(n_steps, B_DIM);
//	TrajArray plan_traj(n_steps, B_DIM+U_DIM);
//	TrajArray exec_open_traj(n_steps, B_DIM);
	double mpc_trans_err, open_trans_err;

	VecBd theta_init = traj.block(0,0,1,B_DIM).transpose();
	VecXd x_init;
	MatXXd rt_Sigma_init;
	brad->decomposeBelief(theta_init, x_init, rt_Sigma_init);
	x_init = brad->mvnrnd(x_init, rt_Sigma_init * rt_Sigma_init.transpose());
	brad->composeBelief(x_init, rt_Sigma_init, theta_init);

	MatrixXd qq(Q_DIM, n_steps);
	for (int j=0; j<n_steps; j++) {
		qq.col(j) = brad->VectorXdRand(Q_DIM);
	}
	MatrixXd rr(R_DIM, n_steps);
	for (int j=0; j<n_steps; j++) {
		rr.col(j) = brad->VectorXdRand(R_DIM);
	}


	VecXd x;
	VecXd x_gt;
	MatXXd rt_Sigma;
	VecBd theta;
	Vector3d trans_gt, trans_est;

	struct timeval startTimeStruct;
	gettimeofday(&startTimeStruct, NULL);
	unsigned long int startTime = startTimeStruct.tv_sec*(long unsigned int)(1e6) + startTimeStruct.tv_usec;

	// MPC
	x = x_init;
	x_gt = x_init;
	rt_Sigma = rt_Sigma_init;
	theta = theta_init;
	int i=0;
	do {
		gPCI = &pci;
		if (root.isMember("costs")) fromJsonArray(root["costs"], pci.cost_infos);
		if (root.isMember("constraints")) fromJsonArray(root["constraints"], pci.cnt_infos);
		gPCI = NULL;

		TrajOptProbPtr prob = ConstructProblem(pci);
		TrajOptResultPtr result = OptimizeProblem(prob, interactive);
		if (i==0) plan_traj = result->traj;
		VecUd u = result->traj.block(0,B_DIM,1,U_DIM).transpose();
		exec_mpc_traj.block(i,0,1,B_DIM) = theta.transpose();
		exec_mpc_traj.block(i,B_DIM,1,U_DIM) = u.transpose();

		// simulate dynamics and observe
		VecQd q = qq.col(i);
		VecRd r = rr.col(i);
		x_gt = brad->Dynamics(x_gt,u,q);
		VecZd z = brad->Observe(x_gt,r);
		exec_mpc_gt_traj.row(i) = x_gt.transpose();

		brad->ekfUpdate(u, x, rt_Sigma, x, rt_Sigma, true, z);

		DblVec lower, upper;
		brad->GetDOFLimits(lower, upper);
		for (int j=0; j<x.size(); j++) {
			if (x(j) < lower[j]) x(j) = lower[j] + STEP;
			if (x(j) > upper[j]) x(j) = upper[j] - STEP;
		}

		brad->composeBelief(x, rt_Sigma, theta);

		brad->SetDOFValues(toDblVec(x));
		cout << "setting robot with dofs values " << endl;
		cout << x.transpose() << endl;
		cout << "actual robot dofs values are " << endl;
		cout << toVectorXd(brad->GetDOFValues()).transpose() << endl;
		traj = result->traj.bottomRows(result->traj.rows()-1);
		traj.block(0,0,1,B_DIM) = theta.transpose();
		pci.basic_info.n_steps--;

		cout << "-------------------------------------------" << endl;

		i++;
	} while (pci.basic_info.n_steps > 1);
	exec_mpc_traj.block(n_steps-1,0,1,B_DIM) = theta.transpose();
	exec_mpc_traj.block(n_steps-1,B_DIM,1,U_DIM) = VecUd::Zero().transpose();
	exec_mpc_gt_traj.row(n_steps-1) = x_gt.transpose();
	mpc_trans_err = (endEffectorPosition(brad, x) - endEffectorPosition(brad, x_gt)).norm();

	gettimeofday(&startTimeStruct, NULL);
	unsigned long int curTime = startTimeStruct.tv_sec*(long unsigned int)(1e6) + startTimeStruct.tv_usec;

	cout << "Total MPC time (s): " << (1e-6) * (curTime - startTime) << endl;

	// non MPC
	x = x_init;
	x_gt = x_init;
	rt_Sigma = rt_Sigma_init;
	theta = theta_init;
	for (int i=0; i<n_steps-1; i++) {
		VecUd u = plan_traj.block(i,B_DIM,1,U_DIM).transpose();

		exec_open_traj.block(i,0,1,B_DIM) = theta.transpose();
		exec_open_traj.block(i,B_DIM,1,U_DIM) = u.transpose();

		VecQd q = qq.col(i);
		VecRd r = rr.col(i);
		x_gt = brad->Dynamics(x_gt,u,q);
		VecZd z = brad->Observe(x_gt,r);
		exec_open_gt_traj.row(i) = x_gt.transpose();

		brad->ekfUpdate(u, x, rt_Sigma, x, rt_Sigma, true, z);
		brad->composeBelief(x, rt_Sigma, theta);
	}
	exec_open_traj.block(n_steps-1,0,1,B_DIM) = theta.transpose();
	exec_open_traj.block(n_steps-1,B_DIM,1,U_DIM) = VecUd::Zero().transpose();
	exec_open_gt_traj.row(n_steps-1) = x_gt.transpose();
	open_trans_err = (endEffectorPosition(brad, x) - endEffectorPosition(brad, x_gt)).norm();

	cout << "-------------------------------------------" << endl;
	cout << "cost exec_mpc_traj " << calcTrajectoryCost(brad, exec_mpc_traj) << "\t" << mpc_trans_err << endl;
	cout << "cost exec_open_traj " << calcTrajectoryCost(brad, exec_open_traj) << "\t" << open_trans_err << endl;
	cout << "cost plan_traj " << calcTrajectoryCost(brad, plan_traj) << endl;

	boost::shared_ptr<CollisionChecker> cc = CollisionChecker::GetOrCreate(*brad->GetRobot()->GetEnv());
	cout << isTrajectoryInCollision(cc, exec_mpc_gt_traj.leftCols(X_DIM), brad) << endl;
	cout << isTrajectoryInCollision(cc, exec_open_gt_traj.leftCols(X_DIM), brad) << endl;
	cout << isTrajectoryInCollision(cc, plan_traj.leftCols(X_DIM), brad) << endl;

	VectorXd stats(8);
	stats << calcTrajectoryCost(brad, exec_mpc_traj), calcTrajectoryCost(brad, exec_open_traj), calcTrajectoryCost(brad, plan_traj), mpc_trans_err, open_trans_err,
			isTrajectoryInCollision(cc, exec_mpc_gt_traj.leftCols(X_DIM), brad), isTrajectoryInCollision(cc, exec_open_gt_traj.leftCols(X_DIM), brad),
			isTrajectoryInCollision(cc, plan_traj.leftCols(X_DIM), brad);

	OSGViewerPtr viewer = OSGViewer::GetOrCreate(brad->GetRobot()->GetEnv());
	vector<GraphHandlePtr> handles;

	bool inter = false;
	if (inter) {
		bool gen_figs = false;
		while (true) {
			cout << "executed MPC" << endl;
			if (gen_figs) dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));
			renderTrajectory(brad, viewer, exec_mpc_traj, handles, 3);
			if (gen_figs) {
				viewer->Idle();
				dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(4,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,0,1));
				viewer->Idle();
			}
			handles.clear();

			cout << "executed open loop" << endl;
			if (gen_figs) dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));
			renderTrajectory(brad, viewer, exec_open_traj, handles, 3);
			if (gen_figs) {
				viewer->Idle();
				dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(4,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,0,1));
				viewer->Idle();
			}
			handles.clear();

			cout << "planned" << endl;
			if (gen_figs) dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));
			renderTrajectory(brad, viewer, plan_traj, handles, 3);
			if (gen_figs) {
				viewer->Idle();
				dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(4,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,0,1));
				viewer->Idle();
			}
			handles.clear();

			cout << "RRT" << endl;
			if (gen_figs) dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(0,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));
			renderTrajectory(brad, viewer, init_traj, handles, 3);
			if (gen_figs) {
				viewer->Idle();
				dynamic_cast<OSGViewer::EventHandler*>(viewer->GetViewer().getCameraManipulator())->setTransformation(osg::Vec3d(4,0,4), osg::Vec3d(0,0,0), osg::Vec3d(0,0,1));
				viewer->Idle();
			}
			handles.clear();
		}
	}
	return stats;
}

TrajOptResultPtr OptimizeProblem(TrajOptProbPtr prob, bool plot) {
	RobotBase::RobotStateSaver saver = prob->GetRAD()->Save();
	BasicTrustRegionSQP opt(prob);
	opt.max_iter_ = 100;
	opt.min_approx_improve_frac_ = .001;
	opt.merit_error_coeff_ = 20;
	opt.max_merit_coeff_increases_ = 10;

	if (plot) opt.addCallback(PlotCallback(*prob));
	//  opt.addCallback(boost::bind(&PlotCosts, boost::ref(prob->getCosts()),boost::ref(*prob->GetRAD()), boost::ref(prob->GetVars()), _1));
	opt.initialize(trajToDblVec(prob->GetInitTraj()));

	struct timeval startTimeStruct;
	gettimeofday(&startTimeStruct, NULL);
	unsigned long int startTime = startTimeStruct.tv_sec*(long unsigned int)(1e6) + startTimeStruct.tv_usec;

	opt.optimize();

	gettimeofday(&startTimeStruct, NULL);
	unsigned long int curTime = startTimeStruct.tv_sec*(long unsigned int)(1e6) + startTimeStruct.tv_usec;

	cout << "Total optimization time (s): " << (1e-6) * (curTime - startTime) << endl;
	TrajOptResultPtr result(new TrajOptResult(opt.results(), *prob));

	return result;
}

TrajOptProbPtr ConstructProblem(const ProblemConstructionInfo& pci) {
	TrajOptProbPtr prob(new TrajOptProb());
	const BasicInfo& bi = pci.basic_info;
	prob->belief_space = bi.belief_space;
	int n_steps = bi.n_steps;

	prob->m_rad = pci.rad;
	int n_dof = prob->m_rad->GetDOF();


	DblVec lower, upper;
	prob->m_rad->GetDOFLimits(lower, upper);
	vector<double> vlower, vupper;
	vector<string> names;
	for (int i=0; i < n_steps; ++i) {
		vlower.insert(vlower.end(), lower.data(), lower.data()+lower.size());
		vupper.insert(vupper.end(), upper.data(), upper.data()+upper.size());
		for (unsigned j=0; j < n_dof; ++j) {
			names.push_back( (boost::format("j_%i_%i")%i%j).str() );
		}
		if (bi.belief_space) {
			for (unsigned jj=0; jj< n_dof; ++jj) {
				for (unsigned ii=jj; ii < n_dof; ++ii) {
					names.push_back( (boost::format("cov_%i_%i_%i")%i%ii%jj).str() );
					vlower.push_back(-INFINITY);
					vupper.push_back(INFINITY);
				}
			}
			for (unsigned j=0; j < n_dof; ++j) {
				names.push_back( (boost::format("u_%i_%i")%i%j).str() );
				vlower.push_back(-INFINITY);
				vupper.push_back(INFINITY);
			}
		}
	}
	prob->createVariables(names, vlower, vupper);

	if (bi.belief_space)
		prob->m_traj_vars = VarArray(n_steps, B_DIM + U_DIM, prob->vars_.data());
	else
		prob->m_traj_vars = VarArray(n_steps, n_dof, prob->vars_.data());

	DblVec cur_dofvals = prob->m_rad->GetDOFValues();

	if (bi.start_fixed) {
		if (pci.init_info.data.rows() > 0 && !allClose(toVectorXd(cur_dofvals), pci.init_info.data.block(0,0,1,n_dof).transpose())) {
			cout << toVectorXd(cur_dofvals).transpose() << endl;
			cout << pci.init_info.data.block(0,0,1,n_dof) << endl;
			LOG_WARN("robot dof values don't match initialization. I don't know what you want me to use for the dof values");
		}
		int n_fixed_terms = n_dof;
		if (bi.belief_space) n_fixed_terms = B_DIM;
		for (int j=0; j < n_fixed_terms; ++j) {
			prob->addLinearConstr(exprSub(AffExpr(prob->m_traj_vars(0,j)), pci.init_info.data(0,j)), EQ);
		}
	}

	if (!bi.dofs_fixed.empty()) {
		BOOST_FOREACH(const int& dof_ind, bi.dofs_fixed) {
			for (int i=1; i < prob->GetNumSteps(); ++i) {
				prob->addLinearConstr(exprSub(AffExpr(prob->m_traj_vars(i,dof_ind)), AffExpr(prob->m_traj_vars(0,dof_ind))), EQ);
			}
		}
	}

	BOOST_FOREACH(const CostInfoPtr& ci, pci.cost_infos) {
		ci->hatch(*prob);
	}
	BOOST_FOREACH(const CntInfoPtr& ci, pci.cnt_infos) {
		ci->hatch(*prob);
	}

	if (bi.belief_space) {
		for (int i=0; i < n_steps-1; i++) {
			VarVector theta0_vars = prob->m_traj_vars.block(0,0,n_steps,B_DIM).row(i);
			VarVector theta1_vars = prob->m_traj_vars.block(0,0,n_steps,B_DIM).row(i+1);
			VarVector u_vars = prob->m_traj_vars.block(0,B_DIM,n_steps,U_DIM).row(i);
			prob->addConstr(ConstraintPtr(new BeliefDynamicsConstraint(theta0_vars, theta1_vars, u_vars, prob->GetRAD())));
		}
	}

	prob->SetInitTraj(pci.init_info.data);

	return prob;

}
TrajOptProbPtr ConstructProblem(const Json::Value& root, OpenRAVE::EnvironmentBasePtr env) {
	ProblemConstructionInfo pci(env);
	pci.fromJson(root);
	return ConstructProblem(pci);
}


TrajOptProb::TrajOptProb(int n_steps, BeliefRobotAndDOFPtr rad) : m_rad(rad) {
	DblVec lower, upper;
	m_rad->GetDOFLimits(lower, upper);
	int n_dof = m_rad->GetDOF();
	vector<double> vlower, vupper;
	vector<string> names;
	for (int i=0; i < n_steps; ++i) {
		vlower.insert(vlower.end(), lower.data(), lower.data()+lower.size());
		vupper.insert(vupper.end(), upper.data(), upper.data()+upper.size());
		for (unsigned j=0; j < n_dof; ++j) {
			names.push_back( (boost::format("j_%i_%i")%i%j).str() );
		}
	}
	createVariables(names, vlower, vupper);
	m_traj_vars = VarArray(n_steps, n_dof, getVars().data());

}


TrajOptProb::TrajOptProb() {
}

void PoseCostInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];
	childFromJson(params, timestep, "timestep", gPCI->basic_info.n_steps-1);
	childFromJson(params, xyz,"xyz");
	childFromJson(params, wxyz,"wxyz");
	childFromJson(params, pos_coeffs,"pos_coeffs", (Vector3d)Vector3d::Ones());
	childFromJson(params, rot_coeffs,"rot_coeffs", (Vector3d)Vector3d::Ones());

	string linkstr;
	childFromJson(params, linkstr, "link");
	link = gPCI->rad->GetRobot()->GetLink(linkstr);
	if (!link) {
		PRINT_AND_THROW(boost::format("invalid link name: %s")%linkstr);
	}
}
CostInfoPtr PoseCostInfo::create() {
	return CostInfoPtr(new PoseCostInfo());
}
void PoseCostInfo::hatch(TrajOptProb& prob) {
	prob.addCost(CostPtr(new CartPoseCost(prob.GetVarRow(timestep), toRaveTransform(wxyz, xyz), rot_coeffs, pos_coeffs, prob.GetRAD(), link)));
	prob.getCosts().back()->setName(name);
}

void PoseCntInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];
	childFromJson(params, timestep, "timestep", gPCI->basic_info.n_steps-1);
	childFromJson(params, xyz,"xyz");
	childFromJson(params, wxyz,"wxyz");
	childFromJson(params, pos_coeffs,"pos_coeffs", (Vector3d)Vector3d::Ones());
	childFromJson(params, rot_coeffs,"rot_coeffs", (Vector3d)Vector3d::Ones());

	string linkstr;
	childFromJson(params, linkstr, "link");
	link = gPCI->rad->GetRobot()->GetLink(linkstr);
	if (!link) {
		PRINT_AND_THROW(boost::format("invalid link name: %s")%linkstr);
	}
}

void JointPosCostInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	int n_steps = gPCI->basic_info.n_steps;
	const Value& params = v["params"];
	childFromJson(params, vals, "vals");
	childFromJson(params, coeffs, "coeffs");
	if (coeffs.size() == 1) coeffs = DblVec(n_steps, coeffs[0]);

	int n_dof = gPCI->rad->GetDOF();
	if (vals.size() != n_dof) {
		PRINT_AND_THROW( boost::format("wrong number of dof vals. expected %i got %i")%n_dof%vals.size());
	}
	childFromJson(params, timestep, "timestep", gPCI->basic_info.n_steps-1);
}
void JointPosCostInfo::hatch(TrajOptProb& prob) {
	prob.addCost(CostPtr(new JointPosCost(prob.GetVarRow(timestep), toVectorXd(vals), toVectorXd(coeffs))));
	prob.getCosts().back()->setName(name);
}
CostInfoPtr JointPosCostInfo::create() {
	return CostInfoPtr(new JointPosCostInfo());
}


CntInfoPtr PoseCntInfo::create() {
	return CntInfoPtr(new PoseCntInfo());
}
void PoseCntInfo::hatch(TrajOptProb& prob) {
	VectorXd coeffs(6); coeffs << rot_coeffs, pos_coeffs;
	prob.addConstr(ConstraintPtr(new CartPoseConstraint(prob.GetVarRow(timestep), toRaveTransform(wxyz, xyz), prob.GetRAD(), link, coeffs)));
	prob.getEqConstraints().back()->setName(name);
}

void CartVelCntInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];
	childFromJson(params, first_step, "first_step");
	childFromJson(params, last_step, "last_step");
	childFromJson(params, distance_limit,"distance_limit");

	FAIL_IF_FALSE((first_step >= 0) && (first_step <= gPCI->basic_info.n_steps-1) && (first_step < last_step));
	FAIL_IF_FALSE((last_step > 0) && (last_step <= gPCI->basic_info.n_steps-1));

	string linkstr;
	childFromJson(params, linkstr, "link");
	link = gPCI->rad->GetRobot()->GetLink(linkstr);
	if (!link) {
		PRINT_AND_THROW( boost::format("invalid link name: %s")%linkstr);
	}
}
CntInfoPtr CartVelCntInfo::create() {
	return CntInfoPtr(new CartVelCntInfo());
}
void CartVelCntInfo::hatch(TrajOptProb& prob) {
	for (int iStep = first_step; iStep < last_step; ++iStep) {
		prob.addConstr(ConstraintPtr(new CartVelConstraint(prob.GetVarRow(iStep), prob.GetVarRow(iStep+1), prob.GetRAD(), link, distance_limit)));
		prob.getIneqConstraints().back()->setName(name);
	}
}

void JointVelCostInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	childFromJson(params, coeffs,"coeffs");
	int n_dof = gPCI->rad->GetDOF();
	if (coeffs.size() == 1) coeffs = DblVec(n_dof, coeffs[0]);
	else if (coeffs.size() != n_dof) {
		PRINT_AND_THROW( boost::format("wrong number of coeffs. expected %i got %i")%n_dof%coeffs.size());
	}
}
CostInfoPtr JointVelCostInfo::create() {
	return CostInfoPtr(new JointVelCostInfo());
}
void JointVelCostInfo::hatch(TrajOptProb& prob) {
	// belief-alex take the submatrix because we want joint-vel only on the joint variables
	prob.addCost(CostPtr(new JointVelCost(prob.GetVars().block(0,0,prob.GetVars().m_nRow, prob.GetRAD()->GetDOF()), toVectorXd(coeffs))));
	prob.getCosts().back()->setName(name);
}

void CollisionCostInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	int n_steps = gPCI->basic_info.n_steps;
	childFromJson(params, coeffs,"coeffs");
	if (coeffs.size() == 1) coeffs = DblVec(n_steps, coeffs[0]);
	else if (coeffs.size() != n_steps) {
		PRINT_AND_THROW( boost::format("wrong size: coeffs. expected %i got %i")%n_steps%coeffs.size() );
	}
	childFromJson(params, dist_pen,"dist_pen");
	if (dist_pen.size() == 1) dist_pen = DblVec(n_steps, dist_pen[0]);
	else if (dist_pen.size() != n_steps) {
		PRINT_AND_THROW( boost::format("wrong size: dist_pen. expected %i got %i")%n_steps%dist_pen.size() );
	}
	childFromJson(params, belief_space,"belief_space", false);
}
void CollisionCostInfo::hatch(TrajOptProb& prob) {
	for (int i=0; i < prob.GetNumSteps(); ++i) {
		if (belief_space)
			prob.addCost(CostPtr(new CollisionCost(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,B_DIM))));
		else
			prob.addCost(CostPtr(new CollisionCost(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,prob.GetRAD()->GetDOF()))));
		prob.getCosts().back()->setName( (boost::format("%s_%i")%name%i).str() );
	}
	CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*prob.GetEnv());
	cc->SetContactDistance(*std::max_element(dist_pen.begin(), dist_pen.end()) + .04);
}
CostInfoPtr CollisionCostInfo::create() {
	return CostInfoPtr(new CollisionCostInfo());
}


void ContinuousCollisionCostInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	int n_steps = gPCI->basic_info.n_steps;
	childFromJson(params, first_step, "first_step", 0);
	childFromJson(params, last_step, "last_step", n_steps-1);
	childFromJson(params, coeffs, "coeffs");
	int n_terms = last_step - first_step;
	if (coeffs.size() == 1) coeffs = DblVec(n_terms, coeffs[0]);
	else if (coeffs.size() != n_terms) {
		PRINT_AND_THROW (boost::format("wrong size: coeffs. expected %i got %i")%n_terms%coeffs.size());
	}
	childFromJson(params, dist_pen,"dist_pen");
	if (dist_pen.size() == 1) dist_pen = DblVec(n_terms, dist_pen[0]);
	else if (dist_pen.size() != n_terms) {
		PRINT_AND_THROW(boost::format("wrong size: dist_pen. expected %i got %i")%n_terms%dist_pen.size());
	}
}
void ContinuousCollisionCostInfo::hatch(TrajOptProb& prob) {
	for (int i=first_step; i < last_step; ++i) {
		prob.addCost(CostPtr(new CollisionCost(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,prob.GetRAD()->GetDOF()), prob.GetVars().rblock(i+1,0,prob.GetRAD()->GetDOF()))));
		prob.getCosts().back()->setName( (boost::format("%s_%i")%name%i).str() );
	}
	CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*prob.GetEnv());
	cc->SetContactDistance(*std::max_element(dist_pen.begin(), dist_pen.end()) + .04);
}
CostInfoPtr ContinuousCollisionCostInfo::create() {
	return CostInfoPtr(new ContinuousCollisionCostInfo());
}


void CollisionCntInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	int n_steps = gPCI->basic_info.n_steps;
	childFromJson(params, coeffs,"coeffs");
	if (coeffs.size() == 1) coeffs = DblVec(n_steps, coeffs[0]);
	else if (coeffs.size() != n_steps) {
		PRINT_AND_THROW( boost::format("wrong size: coeffs. expected %i got %i")%n_steps%coeffs.size() );
	}
	childFromJson(params, dist_pen,"dist_pen");
	if (dist_pen.size() == 1) dist_pen = DblVec(n_steps, dist_pen[0]);
	else if (dist_pen.size() != n_steps) {
		PRINT_AND_THROW( boost::format("wrong size: dist_pen. expected %i got %i")%n_steps%dist_pen.size() );
	}
	childFromJson(params, belief_space,"belief_space", false);
}
void CollisionCntInfo::hatch(TrajOptProb& prob) {
	for (int i=0; i < prob.GetNumSteps(); ++i) {
		if (belief_space)
			prob.addConstr(ConstraintPtr(new CollisionConstraint(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,B_DIM))));
		else
			prob.addConstr(ConstraintPtr(new CollisionConstraint(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,prob.GetRAD()->GetDOF()))));
		prob.getConstraints().back()->setName( (boost::format("%s_%i")%name%i).str() );
	}
	CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*prob.GetEnv());
	cc->SetContactDistance(*std::max_element(dist_pen.begin(), dist_pen.end()) + .04);
}
CntInfoPtr CollisionCntInfo::create() {
	return CntInfoPtr(new CollisionCntInfo());
}


void ContinuousCollisionCntInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	int n_steps = gPCI->basic_info.n_steps;
	childFromJson(params, first_step, "first_step", 0);
	childFromJson(params, last_step, "last_step", n_steps-1);
	childFromJson(params, coeffs, "coeffs");
	int n_terms = last_step - first_step;
	if (coeffs.size() == 1) coeffs = DblVec(n_terms, coeffs[0]);
	else if (coeffs.size() != n_terms) {
		PRINT_AND_THROW (boost::format("wrong size: coeffs. expected %i got %i")%n_terms%coeffs.size());
	}
	childFromJson(params, dist_pen,"dist_pen");
	if (dist_pen.size() == 1) dist_pen = DblVec(n_terms, dist_pen[0]);
	else if (dist_pen.size() != n_terms) {
		PRINT_AND_THROW(boost::format("wrong size: dist_pen. expected %i got %i")%n_terms%dist_pen.size());
	}
}
void ContinuousCollisionCntInfo::hatch(TrajOptProb& prob) {
	for (int i=first_step; i < last_step; ++i) {
		prob.addConstr(ConstraintPtr(new CollisionConstraint(dist_pen[i], coeffs[i], prob.GetRAD(), prob.GetVars().rblock(i,0,prob.GetRAD()->GetDOF()), prob.GetVars().rblock(i+1,0,prob.GetRAD()->GetDOF()))));
		prob.getConstraints().back()->setName( (boost::format("%s_%i")%name%i).str() );
	}
	CollisionCheckerPtr cc = CollisionChecker::GetOrCreate(*prob.GetEnv());
	cc->SetContactDistance(*std::max_element(dist_pen.begin(), dist_pen.end()) + .04);
}
CntInfoPtr ContinuousCollisionCntInfo::create() {
	return CntInfoPtr(new ContinuousCollisionCntInfo());
}


void JointConstraintInfo::fromJson(const Value& v) {
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];
	childFromJson(params, vals, "vals");

	int n_dof = gPCI->rad->GetDOF();
	if (vals.size() != n_dof) {
		PRINT_AND_THROW( boost::format("wrong number of dof vals. expected %i got %i")%n_dof%vals.size());
	}
	childFromJson(params, timestep, "timestep", gPCI->basic_info.n_steps-1);
}

void JointConstraintInfo::hatch(TrajOptProb& prob) {
	VarVector vars = prob.GetVarRow(timestep);
	int n_dof = vars.size();
	for (int j=0; j < n_dof; ++j) {
		prob.addLinearConstr(exprSub(AffExpr(vars[j]), vals[j]), EQ);
	}
}
CntInfoPtr JointConstraintInfo::create() {
	return CntInfoPtr(new JointConstraintInfo());
}



void ControlCostInfo::fromJson(const Value& v) {
	belief_space = gPCI->basic_info.belief_space;
	if (!belief_space) {
		LOG_WARN("control cost can only be used in belief space. ignoring.");
		return;
	}
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	childFromJson(params, coeffs,"coeffs");
	int n_dof = gPCI->rad->GetDOF();
	if (coeffs.size() == 1) coeffs = DblVec(n_dof, coeffs[0]);
	else if (coeffs.size() != n_dof) {
		PRINT_AND_THROW( boost::format("wrong number of coeffs. expected %i got %i")%n_dof%coeffs.size());
	}
}
CostInfoPtr ControlCostInfo::create() {
	return CostInfoPtr(new ControlCostInfo());
}
void ControlCostInfo::hatch(TrajOptProb& prob) {
	if (!belief_space) return;
	prob.addCost(CostPtr(new ControlCost(prob.GetVars().block(0,B_DIM,prob.GetVars().m_nRow-1, prob.GetRAD()->GetDOF()), toVectorXd(coeffs))));
}

void ControlCntInfo::fromJson(const Value& v) {
	belief_space = gPCI->basic_info.belief_space;
	if (!belief_space) {
		LOG_WARN("control constraint can only be used in belief space. ignoring.");
		return;
	}
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];
	childFromJson(params, u_min, "u_min");
	childFromJson(params, u_max, "u_max");
}
CntInfoPtr ControlCntInfo::create() {
	return CntInfoPtr(new ControlCntInfo());
}
void ControlCntInfo::hatch(TrajOptProb& prob) {
	if (!belief_space) return;
	for (int i=0; i < prob.GetVars().m_nRow-1; i++) {
		for (int j=B_DIM; j < B_DIM+U_DIM; j++) {
			prob.addLinearConstr(exprSub(AffExpr(prob.GetVars()(i,j)), u_max), INEQ);
			prob.addLinearConstr(exprMult(exprSub(AffExpr(prob.GetVars()(i,j)), u_min), -1), INEQ);
		}
	}
}

void CovarianceCostInfo::fromJson(const Value& v) {
	belief_space = gPCI->basic_info.belief_space;
	if (!belief_space) {
		LOG_WARN("covariance cost can only be used in belief space. ignoring.");
		return;
	}
	FAIL_IF_FALSE(v.isMember("params"));
	const Value& params = v["params"];

	FAIL_IF_FALSE(params.isMember("Q"));
	const Value& Q_array = params["Q"];
	Json::fromJson(Q_array, Q);

	int n_dof = gPCI->rad->GetDOF();
	if (Q.rows()!=n_dof) PRINT_AND_THROW("cost matrix for the covariance has wrong number of rows");
	if (Q.cols()!=n_dof) PRINT_AND_THROW("cost matrix for the covariance has wrong number of cols");
}
CostInfoPtr CovarianceCostInfo::create() {
	return CostInfoPtr(new CovarianceCostInfo());
}
void CovarianceCostInfo::hatch(TrajOptProb& prob) {
	if (!belief_space) return;
	int n_steps = prob.GetVars().m_nRow;
	for (int i=0; i < n_steps; ++i) {
		VarVector rtSigma_vars = prob.GetVars().rblock(i,U_DIM,B_DIM-U_DIM);
		prob.addCost(CostPtr(new CovarianceCost(rtSigma_vars, Q, prob.GetRAD())));
	}
}

}
