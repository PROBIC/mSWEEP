// Riemannian conjugate gradient for parameter estimation.
#include "rcg.hpp"

#include <assert.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <iostream>

double digamma(double x) {
  double result = 0, xx, xx2, xx4;
  assert(x > 0);
  for ( ; x < 7; ++x)
    result -= 1/x;
  x -= 1.0/2.0;
  xx = 1.0/x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
  return result;
}

void logsumexp(Matrix<double> &gamma_Z, Matrix<double> &q_Z) {
  unsigned n_cols = gamma_Z.get_cols();
  unsigned n_rows = gamma_Z.get_rows();
  for (unsigned i = 0; i < n_cols; ++i) {
    double m = gamma_Z.log_sum_exp_col(i);
    for (unsigned j = 0; j < n_rows; ++j) {
      gamma_Z(j, i) -= m;
      q_Z(j, i) = std::exp(gamma_Z(j, i));
    }
  }
}

double mixt_negnatgrad(const Matrix<double> &q_Z, const Matrix<double> &gamma_Z, const std::vector<double> &N_k, const Matrix<double> &logl, Matrix<double> &dL_dphi) {
  std::vector<double> colsums(q_Z.get_cols());  
  for (unsigned i = 0; i < dL_dphi.get_rows(); ++i) {
    double digamma_N_k = digamma(N_k[i]);
    for (unsigned j = 0; j < dL_dphi.get_cols(); ++j) {
      dL_dphi(i, j) = logl(i, j);
      dL_dphi(i, j) += digamma_N_k - gamma_Z(i, j) - 1;
      colsums[j] += dL_dphi(i, j) * q_Z(i, j);
    }
  }
  double newnorm = 0.0;
  for (unsigned i = 0; i < dL_dphi.get_rows(); ++i) {
    for (unsigned j = 0; j < dL_dphi.get_cols(); ++j) {
      // dL_dgamma(i, j) would be q_Z(i, j) * (dL_dphi(i, j) - colsums[j])
      newnorm += q_Z(i, j) * (dL_dphi(i, j) - colsums[j]) * dL_dphi(i, j);
    }
  }
  return newnorm;
}

void ELBO_rcg_mat(const Matrix<double> &logl, const Matrix<double> &q_Z, const Matrix<double> &gamma_Z, const std::vector<long unsigned> &counts, const std::vector<double> &alpha0, const std::vector<double> &N_k, long double &bound) {
  for (unsigned i = 0; i < q_Z.get_rows(); ++i) {
    for (unsigned j = 0; j < q_Z.get_cols(); ++j) {
      bound += (q_Z(i, j) * (logl(i, j) - gamma_Z(i, j))) * counts[j];
    }
    bound -= lgamma(alpha0[i]) - lgamma(N_k[i]);
  }
}

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const Sample &sample, const std::vector<double> &alpha0, const double &tol, unsigned maxiters) {
  unsigned n_rows = logl.get_rows();
  unsigned n_cols = logl.get_cols();
  Matrix<double> gamma_Z(n_rows, n_cols, 1.0);
  Matrix<double> gamma_new(n_rows, n_cols, 0.0);
  Matrix<double> oldstep(n_rows, n_cols, 0.0);
  Matrix<double> step(n_rows, n_cols, 0.0);
  double oldnorm = 1.0;

  Matrix<double> q_Z(n_rows, n_cols, 0.0);
  logsumexp(gamma_Z, q_Z);
  long double bound = -100000.0;
  bool didreset = false;

  double bound_const = std::accumulate(alpha0.begin(), alpha0.end(), 0.0);;
  bound_const += sample.total_counts();
  bound_const = -lgamma(bound_const);
  bound_const += std::accumulate(alpha0.begin(), alpha0.end(), 0.0, [](double a, double b) { return a + lgamma(b); });

  std::vector<double> N_k(alpha0.size());
  q_Z.right_multiply(sample.ec_counts, N_k);
  std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());

  for (unsigned k = 0; k < maxiters; ++k) {
    double newnorm = mixt_negnatgrad(q_Z, gamma_Z, N_k, logl, step);
    double beta_FR = newnorm/oldnorm;
    oldnorm = newnorm;

    if (didreset) {
      oldstep *= 0.0;
    } else if (beta_FR > 0) {
      oldstep *= beta_FR;
      step += oldstep;
    }
    didreset = false;
    
    gamma_new.sum_fill(gamma_Z, step); // gamma_new = gamma_Z + step
    logsumexp(gamma_new, q_Z);
    q_Z.right_multiply(sample.ec_counts, N_k);
    std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());
    long double oldbound = bound;
    bound = bound_const;
    ELBO_rcg_mat(logl, q_Z, gamma_new, sample.ec_counts, alpha0, N_k, bound);

    if (bound < oldbound) {
      didreset = true;
      gamma_Z += step;
      if (beta_FR > 0) {
	gamma_Z -= oldstep;
      }
      logsumexp(gamma_Z, q_Z);
      q_Z.right_multiply(sample.ec_counts, N_k);
      std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());
      bound = bound_const;
      ELBO_rcg_mat(logl, q_Z, gamma_Z, sample.ec_counts, alpha0, N_k, bound);
    } else {
      oldstep = step;
      gamma_Z = gamma_new;
    }
    if (k % 5 == 0) {
      std::cerr << "  " <<  "iter: " << k << ", bound: " << bound << ", |g|: " << newnorm << std::endl;
    }
    if (bound - oldbound < tol && !didreset) {
      logsumexp(gamma_Z, gamma_Z);
      return(gamma_Z);
    }
  }
  logsumexp(gamma_Z, gamma_Z);
  return(gamma_Z);
}

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const long unsigned &total_counts, const std::vector<long unsigned> &ec_counts, const std::vector<double> &alpha0, const double &tol, unsigned maxiters) {
  unsigned n_rows = logl.get_rows();
  unsigned n_cols = logl.get_cols();
  Matrix<double> gamma_Z(n_rows, n_cols, 1.0);
  Matrix<double> gamma_new(n_rows, n_cols, 0.0);
  Matrix<double> oldstep(n_rows, n_cols, 0.0);
  Matrix<double> step(n_rows, n_cols, 0.0);
  double oldnorm = 1.0;

  Matrix<double> q_Z(n_rows, n_cols, 0.0);
  logsumexp(gamma_Z, q_Z);
  long double bound = -100000.0;
  bool didreset = false;

  double bound_const = std::accumulate(alpha0.begin(), alpha0.end(), 0.0);;
  bound_const += total_counts;
  bound_const = -lgamma(bound_const);
  bound_const += std::accumulate(alpha0.begin(), alpha0.end(), 0.0, [](double a, double b) { return a + lgamma(b); });

  std::vector<double> N_k(alpha0.size());
  q_Z.right_multiply(ec_counts, N_k);
  std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());

  for (unsigned k = 0; k < maxiters; ++k) {
    double newnorm = mixt_negnatgrad(q_Z, gamma_Z, N_k, logl, step);
    double beta_FR = newnorm/oldnorm;
    oldnorm = newnorm;

    if (didreset) {
      oldstep *= 0.0;
    } else if (beta_FR > 0) {
      oldstep *= beta_FR;
      step += oldstep;
    }
    didreset = false;
    
    gamma_new.sum_fill(gamma_Z, step); // gamma_new = gamma_Z + step
    logsumexp(gamma_new, q_Z);
    q_Z.right_multiply(ec_counts, N_k);
    std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());
    long double oldbound = bound;
    bound = bound_const;
    ELBO_rcg_mat(logl, q_Z, gamma_new, ec_counts, alpha0, N_k, bound);

    if (bound < oldbound) {
      didreset = true;
      gamma_Z += step;
      if (beta_FR > 0) {
	gamma_Z -= oldstep;
      }
      logsumexp(gamma_Z, q_Z);
      q_Z.right_multiply(ec_counts, N_k);
      std::transform(N_k.begin(), N_k.end(), alpha0.begin(), N_k.begin(), std::plus<double>());
      bound = bound_const;
      ELBO_rcg_mat(logl, q_Z, gamma_Z, ec_counts, alpha0, N_k, bound);
    } else {
      oldstep = step;
      gamma_Z = gamma_new;
    }
    if (k % 5 == 0) {
      std::cerr << "  " <<  "iter: " << k << ", bound: " << bound << ", |g|: " << newnorm << std::endl;
    }
    if (bound - oldbound < tol && !didreset) {
      logsumexp(gamma_Z, gamma_Z);
      return(gamma_Z);
    }
  }
  logsumexp(gamma_Z, gamma_Z);
  return(gamma_Z);
}
