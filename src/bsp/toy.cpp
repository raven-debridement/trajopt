#include "bsp.hpp"

using namespace BSP;

template<typename T, size_t N>
T* end(T (&ra)[N]) {
  return ra + N;
}

class ToyStateFunc : public StateFunc {
public:
  ToyStateFunc(BSPProblemHelperPtr helper);
  virtual int output_size() const;
  virtual VectorXd operator()(const VectorXd& x, const VectorXd& u, const VectorXd& m) const ;
};

class ToyObserveFunc : public ObserveFunc {
public:
  ToyObserveFunc(BSPProblemHelperPtr helper);
  virtual int output_size() const;
  virtual VectorXd operator()(const VectorXd& x, const VectorXd& n) const;
};

class ToyBeliefFunc : public BeliefFunc {
public:
  ToyBeliefFunc(BSPProblemHelperPtr helper);
  ToyBeliefFunc(BSPProblemHelperPtr helper, double alpha);
  virtual int output_size() const;
  virtual VectorXd operator()(const VectorXd& b, const VectorXd& u, const StateFunc& f) const;
  void set_alpha(double alpha);
  void set_tol(double tol);
protected:
  double alpha;
  double tol;
  MatrixXd get_gamma(const VectorXd& x) const;
};

VectorXd ToyStateFunc::operator()(const VectorXd& x // state
                                , const VectorXd& u // control
                                , const VectorXd& m // state noise
                                 ) const {
  assert(x.size() == helper->state_dim);
  assert(u.size() == helper->control_dim);
  assert(m.size() == helper->state_noise_dim);
  VectorXd new_x(helper->state_dim);
  for (int i = 0; i < helper->state_dim; ++i) {
    new_x(i) = x(i) + u(i)*helper->input_dt + sqrt(0.01*u(i)*u(i)*input_dt + 0.001) * m(i);
  }
  return new_x;
}

VectorXd ToyObserveFunc::operator()(const VectorXd& x // state
                                  , const VectorXd& n // observe noise
                                   ) const {
  assert(x.size() == helper->state_dim);
  assert(n.size() == helper->observe_noise_dim);
  return x + 0.1 * n;
}

VectorXd ToyBeliefFunc::operator()(const VectorXd& b    // current belief
                                 , const VectorXd& u    // control
                                 , const StateFunc& f   // state propagation function xt -> xt+1
                                 , const ObserveFunc& h // observation function xt -> yt
                                  ) const {

  assert (b.size() == helper->belief_dim);
  assert (u.size() == helper->control_dim);

  VectorXd x, new_x, new_b;
  MatrixXd sigma, A, M, gamma, H, N, K;

  extract_state(b, &x);
  extract_sigma(b, &sigma);

  // linearize the state propagation function around current point
  f.linearize(x, u, VectorXd::Zero(helper->state_noise_dim), &A, NULL, &M);

  sigma = A*sigma*A.transpose() + M*M.transpose();

  new_x = f.call(x, u, VectorXd::Zero(helper->state_noise_dim));

  gamma = get_gamma(new_x); // TODO make this more generalizable

  // linearize the observation function around current point
  h.linearize(new_x, VectorXd::Zero(helper->observe_noise_dim), &H, &N); 

  H = gamma * H;

  PartialPivLU<MatrixXd> solver();
  K = matrix_div(sigma*H.transpose(), H*sigma* H.transpose() + N*N.transpose());//           solver.solve( ().transpose ).transpose();

  sigma = sigma - gamma*K*(H*sigma);

  return compose_belief(new_x, matrix_sqrt(sigma));
}

bool ToyBeliefFunc::sgndist(const Vector2d& x, Vector2d* dists) const {
  Vector2d p1; p1 << 0, 2;
  Vector2d p2; p2 << 0, 0;
  *dists(0) = (x - p1).norm() - 0.5;
  *dists(1) = (x - p2).norm() - 0.5;
  return (*dists(0) < 0) || (*dists(1) < 0);
}

MatrixXd ToyBeliefFunc::get_gamma(const VectorXd& x) const {
  Vector2d dists;
  sgndist(x, dists);
  double gamma1 = 1. - (1./(1.+exp(-alpha*(dists(0)+tol))));
  double gamma2 = 1. - (1./(1.+exp(-alpha*(dists(1)+tol))));
  Matrix2d gamma;
  gamma << gamma1, 0,
           0, gamma2;
  return gamma;
}

int main(int argc, char** argv) {
  int T = 20;
  int n_dof = 2;

  double start_vec_array[] = {-5, 2};
  double goal_vec_array[] = {-5, 1};

  vector<double> start_vec(start_vec_array, end(start_vec_array));
  vector<double> goal_vec(goal_vec_array, end(goal_vec_array));

  {
    Config config;
    config.add(new Parameter< vector<double> >("s", &start_vec, "s"));
    config.add(new Parameter< vector<double> >("g", &goal_vec, "g"));
    CommandParser parser(config);
    parser.read(argc, argv);
  }

  Vector2d start = toVectorXd(start_vec);
  Vector2d goal = toVectorXd(goal_vec);

  //B b;

  //cout << boost::bind(&A::a, b)() << endl;

  //OptProbPtr prob(new OptProb());

  //BSPHelperPtr helper(new BSPHelper());
  //helper->start = start;
  //helper->goal = goal;
  //helper->T = T;
  //helper->n_dof = n_dof;
  //helper->ConfigureProblem(*prob);

  //BSPTrustRegionSQP opt(prob);
  //opt.max_iter_ = 500;    

  //helper->ConfigureOptimizer(*prob, opt);

  //opt.optimize();

  cout << start.transpose() << endl;
  cout << goal.transpose() << endl;
  return 0;  
}
