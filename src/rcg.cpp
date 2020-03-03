// Riemannian conjugate gradient for parameter estimation.
#include "rcg.hpp"
#include "openmp_config.hpp"

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#endif

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#endif

double digamma(double x) {
  double result = 0, xx, xx2, xx4;
  assert(x > 0);
  for ( ; x < 7; ++x)
    result -= 1/x;
  x -= 1.0/2.0;
  xx = 1.0/x;
  xx2 = xx*xx;
  xx4 = xx2*xx2;
  result += std::log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
  return result;
}

void logsumexp(Matrix<double> &gamma_Z) {
  unsigned n_cols = gamma_Z.get_cols();
  unsigned short n_rows = gamma_Z.get_rows();

  std::vector<double> m(n_cols, 0.0);
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < n_cols; ++i) {
    m[i] = gamma_Z.log_sum_exp_col(i);
  }

#pragma omp parallel for schedule(static)
  for (unsigned short i = 0; i < n_rows; ++i) {
    for (unsigned j = 0; j < n_cols; ++j) {
      gamma_Z(i, j) -= m[j];
    }
  }
}

double mixt_negnatgrad(const Matrix<double> &gamma_Z, const std::vector<double> &N_k, const Matrix<double> &logl, Matrix<double> &dL_dphi) {
  unsigned n_cols = gamma_Z.get_cols();
  unsigned short n_rows = gamma_Z.get_rows();

  std::vector<double> colsums(n_cols, 0.0);
#pragma omp parallel for schedule(static) reduction(vec_double_plus:colsums)
  for (unsigned short i = 0; i < n_rows; ++i) {
    double digamma_N_k = digamma(N_k[i]) - 1.0;
    for (unsigned j = 0; j < n_cols; ++j) {
      dL_dphi(i, j) = logl(i, j);
      dL_dphi(i, j) += digamma_N_k - gamma_Z(i, j);
      colsums[j] += dL_dphi(i, j) * std::exp(gamma_Z(i, j));
    }
  }

  double newnorm = 0.0;
#pragma omp parallel for schedule(static) reduction(+:newnorm)
  for (unsigned short i = 0; i < n_rows; ++i) {
    for (unsigned j = 0; j < n_cols; ++j) {
      // dL_dgamma(i, j) would be q_Z(i, j) * (dL_dphi(i, j) - colsums[j])
      // newnorm += q_Z(i, j) * (dL_dphi(i, j) - colsums[j]) * dL_dphi(i, j);
      newnorm += std::exp(gamma_Z(i, j)) * (dL_dphi(i, j) - colsums[j]) * dL_dphi(i, j);
    }
  }
  return newnorm;
}

void ELBO_rcg_mat(const Matrix<double> &logl, const Matrix<double> &gamma_Z, const std::vector<double> &counts, const std::vector<double> &alpha0, const std::vector<double> &N_k, long double &bound) {
  unsigned short n_rows = gamma_Z.get_rows();
  unsigned n_cols = gamma_Z.get_cols();
#pragma omp parallel for schedule(static) reduction(+:bound)
  for (unsigned short i = 0; i < n_rows; ++i) {
    for (unsigned j = 0; j < n_cols; ++j) {
      bound += std::exp(gamma_Z(i, j) + counts[j])*(logl(i, j) - gamma_Z(i, j));
      //      bound += (std::exp(gamma_Z(i, j)) * (logl(i, j) - gamma_Z(i, j))) * counts[j];
    }
    bound -= std::lgamma(alpha0[i]) - std::lgamma(N_k[i]);
  }
}

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const Sample &sample, const std::vector<double> &alpha0, const double &tol, unsigned maxiters) {
  unsigned short n_rows = logl.get_rows();
  unsigned n_cols = logl.get_cols();
  Matrix<double> gamma_Z(n_rows, n_cols, std::log(1.0/(double)n_rows)); // where gamma_Z is init at 1.0
  Matrix<double> gamma_new(n_rows, n_cols, 0.0);
  Matrix<double> oldstep(n_rows, n_cols, 0.0);
  Matrix<double> step(n_rows, n_cols, 0.0);
  double oldnorm = 1.0;
  long double bound = -100000.0;
  bool didreset = false;
  double bound_const = sample.total_counts();

#pragma omp parallel for schedule(static) reduction(+:bound_const)
  for (unsigned short i = 0; i < n_rows; ++i) {
    bound_const += alpha0[i];
    bound_const += std::lgamma(alpha0[i]);
  }

  bound_const = -std::lgamma(bound_const);
  std::vector<double> N_k(alpha0.size());
  gamma_Z.exp_right_multiply(sample.ec_counts, N_k);

#pragma omp parallel for schedule(static)
    for (unsigned short i = 0; i < n_rows; ++i) {
      N_k[i] += alpha0[i];
    }

  for (unsigned short k = 0; k < maxiters; ++k) {
    double newnorm = mixt_negnatgrad(gamma_Z, N_k, logl, step);
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
    logsumexp(gamma_new);
    gamma_new.exp_right_multiply(sample.ec_counts, N_k);

    #pragma omp parallel for schedule(static)
    for (unsigned short i = 0; i < n_rows; ++i) {
      N_k[i] += alpha0[i];
    }

    long double oldbound = bound;
    bound = bound_const;
    ELBO_rcg_mat(logl, gamma_new, sample.ec_counts, alpha0, N_k, bound);

    if (bound < oldbound) {
      didreset = true;
      gamma_Z += step;
      if (beta_FR > 0) {
	gamma_Z -= oldstep;
      }
      logsumexp(gamma_Z);
      gamma_Z.exp_right_multiply(sample.ec_counts, N_k);

    #pragma omp parallel for schedule(static)
    for (unsigned short i = 0; i < n_rows; ++i) {
      N_k[i] += alpha0[i];
    }

      bound = bound_const;
      ELBO_rcg_mat(logl, gamma_Z, sample.ec_counts, alpha0, N_k, bound);
    } else {
      oldstep = step;
      gamma_Z = gamma_new;
    }
    if (k % 5 == 0) {
      std::cerr << "  " <<  "iter: " << k << ", bound: " << bound << ", |g|: " << newnorm << std::endl;
    }
    if (bound - oldbound < tol && !didreset) {
      logsumexp(gamma_Z);
// #pragma omp parallel for schedule(static)
//       for (unsigned short i = 0; i < n_rows; ++i) {
// 	for (unsigned j = 0; j < n_cols; ++j) {
// 	  gamma_Z(i, j) = std::exp(gamma_Z(i, j));
// 	}
//       }
      return(gamma_Z);
    }
  }
  logsumexp(gamma_Z);
// #pragma omp parallel for schedule(static)
//       for (unsigned short i = 0; i < n_rows; ++i) {
// 	for (unsigned j = 0; j < n_cols; ++j) {
// 	  gamma_Z(i, j) = std::exp(gamma_Z(i, j));
// 	}
//       }
  return(gamma_Z);
}

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const long unsigned &total_counts, const std::vector<double> &ec_counts, const std::vector<double> &alpha0, const double &tol, unsigned maxiters) {
  unsigned short n_rows = logl.get_rows();
  unsigned n_cols = logl.get_cols();
  Matrix<double> gamma_Z(n_rows, n_cols, std::log(1.0/(double)n_rows)); // where gamma_Z is init at 1.0
  Matrix<double> gamma_new(n_rows, n_cols, 0.0);
  Matrix<double> oldstep(n_rows, n_cols, 0.0);
  Matrix<double> step(n_rows, n_cols, 0.0);
  double oldnorm = 1.0;
  long double bound = -100000.0;
  bool didreset = false;

  double bound_const = total_counts;

  #pragma omp parallel for schedule(static) reduction(+:bound_const)
  for (unsigned short i = 0; i < n_rows; ++i) {
    bound_const += alpha0[i];
    bound_const += std::lgamma(alpha0[i]);
  }

  bound_const = -std::lgamma(bound_const);
  std::vector<double> N_k(alpha0.size());
  gamma_Z.exp_right_multiply(ec_counts, N_k);

  #pragma omp parallel for schedule(static)
  for (unsigned short i = 0; i < n_rows; ++i) {
    N_k[i] += alpha0[i];
  }

  for (unsigned short k = 0; k < maxiters; ++k) {
    double newnorm = mixt_negnatgrad(gamma_Z, N_k, logl, step);
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
    logsumexp(gamma_new);
    //    q_Z.right_multiply(ec_counts, N_k);
    gamma_new.exp_right_multiply(ec_counts, N_k);

    #pragma omp parallel for schedule(static)
    for (unsigned short i = 0; i < n_rows; ++i) {
      N_k[i] += alpha0[i];
    }

    long double oldbound = bound;
    bound = bound_const;
    ELBO_rcg_mat(logl, gamma_new, ec_counts, alpha0, N_k, bound);

    if (bound < oldbound) {
      didreset = true;
      gamma_Z += step;
      if (beta_FR > 0) {
	gamma_Z -= oldstep;
      }
      logsumexp(gamma_Z);
      //      q_Z.right_multiply(ec_counts, N_k);
      gamma_Z.exp_right_multiply(ec_counts, N_k);

      #pragma omp parallel for schedule(static)
      for (unsigned short i = 0; i < n_rows; ++i) {
        N_k[i] += alpha0[i];
      }

      bound = bound_const;
      ELBO_rcg_mat(logl, gamma_Z, ec_counts, alpha0, N_k, bound);
    } else {
      oldstep = step;
      gamma_Z = gamma_new;
    }
    if (k % 5 == 0) {
      std::cerr << "  " <<  "iter: " << k << ", bound: " << bound << ", |g|: " << newnorm << std::endl;
    }
    if (bound - oldbound < tol && !didreset) {
      logsumexp(gamma_Z);
#pragma omp parallel for schedule(static)
      for (unsigned short i = 0; i < n_rows; ++i) {
	for (unsigned j = 0; j < n_cols; ++j) {
	  gamma_Z(i, j) = std::exp(gamma_Z(i, j));
	}
      }
      return(gamma_Z);
    }
  }
  logsumexp(gamma_Z);
#pragma omp parallel for schedule(static)
      for (unsigned short i = 0; i < n_rows; ++i) {
	for (unsigned j = 0; j < n_cols; ++j) {
	  gamma_Z(i, j) = std::exp(gamma_Z(i, j));
	}
      }
  return(gamma_Z);
}
