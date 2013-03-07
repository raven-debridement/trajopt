#include "expr_ops.hpp"
#include "expr_op_overloads.hpp"
#include <cmath>

static inline double sq(double x) {return x*x;}

namespace ipi {
namespace sco {
QuadExpr exprSquare(const Var& a) {
  QuadExpr out;
  out.coeffs.push_back(1);
  out.vars1.push_back(a);
  out.vars2.push_back(a);
  return out;
}

QuadExpr exprSquare(const AffExpr& affexpr) {
  QuadExpr out;
  size_t naff = affexpr.coeffs.size();
  size_t nquad = (naff*(naff+1))/2;

  out.affexpr.constant = sq(affexpr.constant);

  out.affexpr.vars = affexpr.vars;
  out.affexpr.coeffs.resize(naff);
  for (size_t i=0; i < naff; ++i) out.affexpr.coeffs[i] = 2*affexpr.constant*affexpr.coeffs[i];

  out.coeffs.reserve(nquad);
  out.vars1.reserve(nquad);
  out.vars2.reserve(nquad);
  for (size_t i=0; i < naff; ++i) {
    out.vars1.push_back(affexpr.vars[i]);
    out.vars2.push_back(affexpr.vars[i]);
    out.coeffs.push_back(sq(affexpr.coeffs[i]));
    for (size_t j=i+1; j < naff; ++j) {
      out.vars1.push_back(affexpr.vars[i]);
      out.vars2.push_back(affexpr.vars[j]);
      out.coeffs.push_back(2 * affexpr.coeffs[i] * affexpr.coeffs[j]);
    }
  }
  return out;
}

QuadExpr exprMult(const Var& a, const Var& b) {
	QuadExpr expr;
	expr.coeffs.push_back(1);
	expr.vars1.push_back(a);
	expr.vars2.push_back(b);
	return expr;
}

BasicArray<QuadExpr> exprMult(const BasicArray<Var>& A, const BasicArray<Var>& B) {
	BasicArray<QuadExpr> C(A.rows(), B.cols());
	assert(A.cols() == B.rows());
	int m = A.cols();
	for (int i=0; i<A.rows(); i++) {
		for (int j=0; j<B.cols(); j++) {
			C(i,j) = QuadExpr();
			for (int k=0; k<m; k++)
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
		}
	}
	return C;
}

BasicArray<QuadExpr> exprMult(const Eigen::MatrixXd& A, const BasicArray<QuadExpr>& B) {
	BasicArray<QuadExpr> C(A.rows(), B.cols());
	assert(A.cols() == B.rows());
	int m = A.cols();
	for (int i=0; i<A.rows(); i++) {
		for (int j=0; j<B.cols(); j++) {
			C(i,j) = QuadExpr();
			for (int k=0; k<m; k++)
				C(i,j) = C(i,j) + A(i,k)*B(k,j);
		}
	}
	return C;
}

AffExpr cleanupAff(const AffExpr& a) {
  AffExpr out;
  for (size_t i=0; i < a.size(); ++i) {
    if (fabs(a.coeffs[i]) > 1e-7) {
      out.coeffs.push_back(a.coeffs[i]);
      out.vars.push_back(a.vars[i]);
    }
  }
  out.constant = a.constant;
  return out;
}

QuadExpr cleanupQuad(const QuadExpr& q) {
  QuadExpr out;
  out.affexpr = cleanupAff(q.affexpr);
  for (size_t i=0; i < q.size(); ++i) {
    if (fabs(q.coeffs[i]) > 1e-8) {
      out.coeffs.push_back(q.coeffs[i]);
      out.vars1.push_back(q.vars1[i]);
      out.vars2.push_back(q.vars2[i]);
    }
  }
  return out;
}

///////////////////////////////////////////////////////////////


}
}
