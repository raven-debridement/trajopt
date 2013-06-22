#include "bsp.hpp"

using namespace BSP;

template<typename T, size_t N>
T* end(T (&ra)[N]) {
  return ra + N;
}

class ToyBSPProblemHelper : public BSPProblemHelper {
public:
  int input_dt;
};

typedef boost::shared_ptr<ToyBSPProblemHelper> ToyBSPProblemHelperPtr;

typedef Matrix<double, 5, 1> Vector5d;

class ToyStateFunc : public StateFunc<Vector2d, Vector2d, Vector2d> {
public:
  ToyStateFunc(ToyBSPProblemHelperPtr helper);
  virtual int output_size() const;
  virtual Vector2d operator()(const Vector2d& x, const Vector2d& u, const Vector2d& m) const ;
protected:
  ToyBSPProblemHelperPtr helper;
};

class ToyObserveFunc : public ObserveFunc<Vector2d, Vector2d, Vector2d> {
public:
  ToyObserveFunc(ToyBSPProblemHelperPtr helper);
  virtual int output_size() const;
  virtual Vector2d operator()(const Vector2d& x, const Vector2d& n) const;
protected:
  ToyBSPProblemHelperPtr helper;
};

class ToyBeliefFunc : public BeliefFunc<ToyStateFunc, ToyObserveFunc> {
public:
  ToyBeliefFunc(ToyBSPProblemHelperPtr helper);
  ToyBeliefFunc(ToyBSPProblemHelperPtr helper, double alpha);
  virtual int output_size() const;
  virtual Vector5d operator()(const Vector5d& b, const Vector2d& u) const;
  void set_alpha(double alpha);
  void set_tol(double tol);
protected:
  double alpha;
  double tol;
  ToyBSPProblemHelperPtr helper;
  Matrix2d get_gamma(const Vector2d& x) const;
  bool sgndist(const Vector2d& x, Vector2d* dists) const;
  //void extract_state(const VectorXd& belief, VectorXd* output_state) const;
  //void extract_sqrt_sigma(const VectorXd& belief, MatrixXd* output_sqrt_sigma) const;
  //void extract_sigma(const VectorXd& belief, MatrixXd* output_sigma) const;
  //void compose_belief(const VectorXd& state, const MatrixXd& sqrt_sigma, VectorXd* output_belief) const;
};

Vector2d ToyStateFunc::operator()(const Vector2d& x // state
                                , const Vector2d& u // control
                                , const Vector2d& m // state noise
                                 ) const {
  Vector2d new_x;
  for (int i = 0; i < helper->state_dim; ++i) {
    new_x(i) = x(i) + u(i)*helper->input_dt + sqrt(0.01*u(i)*u(i)*helper->input_dt + 0.001) * m(i);
  }
  return new_x;
}

Vector2d ToyObserveFunc::operator()(const Vector2d& x // state
                                  , const Vector2d& n // observe noise
                                   ) const {
  return x + 0.1 * n;
}



Vector5d ToyBeliefFunc::operator()(const Vector5d& b    // current belief
                                 , const Vector2d& u    // control
                                  ) const {

  Vector2d x, new_x;
  Vector5d new_b;
  Matrix2d sigma, A, M, gamma, H, N, K;

  extract_state(b, &x);
  extract_sigma(b, &sigma);

  // linearize the state propagation function around current point

  Vector2d zero = Vector2d::Zero();

  f->linearize(x, u, zero, &A, (ToyStateFunc::ControlGradT*) NULL, &M);

  sigma = A*sigma*A.transpose() + M*M.transpose();

  new_x = f->call(x, u, zero);

  gamma = get_gamma(new_x); // TODO make this more generalizable

  // linearize the observation function around current point
  h->linearize(new_x, zero, &H, &N); 

  H = gamma * H;

  PartialPivLU<Matrix2d> solver();
  K = matrix_div((Matrix2d) (sigma*H.transpose()), (Matrix2d) (H*sigma* H.transpose() + N*N.transpose()));//           solver.solve( ().transpose ).transpose();

  sigma = sigma - gamma*K*(H*sigma);

  compose_belief(new_x, matrix_sqrt(sigma), &new_b);

  return new_b;
}

bool ToyBeliefFunc::sgndist(const Vector2d& x, Vector2d* dists) const {
  Vector2d p1; p1 << 0, 2;
  Vector2d p2; p2 << 0, 0;
  (*dists)(0) = (x - p1).norm() - 0.5;
  (*dists)(1) = (x - p2).norm() - 0.5;
  return (*dists)(0) < 0 || (*dists)(1) < 0;
}

Matrix2d ToyBeliefFunc::get_gamma(const Vector2d& x) const {
  Vector2d dists;
  sgndist(x, &dists);
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
