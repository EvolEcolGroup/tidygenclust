// The fitting loop of gt_radmix() for one value of K: a starting point, a
// few EM steps to prime it, then block updates accelerated by quasi-Newton
// steps until the log-likelihood stops improving. See
// dev_notes/admixture_implementation.md for the algorithm, and radmix_kernels.h
// for the model and the layout of q (K x I) and p (K x J).
//
// The whole loop runs in C++, with every estimate in a buffer allocated once
// at the start of the fit: accepting an estimate swaps buffers, so a fit
// holds a fixed number of copies of the parameters, 3 + 2 n_qn (x0, x1, x2
// and the secant pairs), whatever the number of iterations (plus the
// matrices returned to R, built at the end). Between passes over the
// genotypes, on the main thread, the loop reports progress (FitLog) and
// checks for user interrupts; nothing inside the passes calls R.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include "radmix_geno_access.h"
#include "radmix_kernels.h"
#include "radmix_quasi_newton.h"

using namespace Rcpp;

namespace radmix {

// x with `digits` significant digits (printf's %g, which drops trailing
// zeros, close to what R's format(x, digits = digits) prints), and NaN, Inf
// and -Inf spelt as in R, for the log
static std::string fmt(double x, int digits) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.*g", digits, x);
  return buf;
}

// The lines of the log, also reported as they come with base::message() if
// verbose (so that suppressMessages() silences them). Only called from the
// main thread, outside parallel regions, where calling R is safe.
struct FitLog {
  // whether to print the lines as they come
  bool verbose;
  // all the lines so far, returned to R as the `log` of the fit
  std::vector<std::string> lines;
  // append a line, and print it if verbose (message() is looked up in the
  // base environment, so that a function of the same name elsewhere cannot
  // be called instead)
  void add(const std::string& line) {
    lines.push_back(line);
    if (verbose) {
      Function message("message", R_BaseEnv);
      message(line);
    }
  }
};

}  // namespace radmix

using namespace radmix;

// Fit the model for one value of K.
//
// data: the genotypes, as prepared by radmix_prepare() (see
//   radmix_geno_access.h), with an optional cross-validation mask.
// k: the number of ancestral populations.
// init_q, init_p: optional starting values in the user layout (individuals
//   x K and loci x K), already validated in R; NULL to start at random.
// fix_q, fix_p: keep the corresponding block at its starting value.
// sup_q: optional I x K matrix for a supervised analysis, already validated
//   in R: the row of each reference individual holds its fixed ancestry
//   proportions (inside the feasible set), and the rows of the individuals
//   to estimate are NA. NULL for an unsupervised analysis.
// max_iter, rtol: maximum number of iterations and convergence threshold on
//   the relative change in log-likelihood.
// n_em_iter: number of EM steps that prime a random start (skipped when both
//   init_q and init_p are given).
// n_qn: number of secant pairs used by the quasi-Newton step.
// n_cores: number of OpenMP threads of the passes over the genotypes (the
//   results do not depend on it).
// verbose: print the lines of the log as they come.
// Returns list(q = I x K, p = J x K, loglik, loglik_history, n_iter,
// converged, log), with the estimates in the user layout; loglik_history
// holds the log-likelihood of the start and of the end of each iteration.
// [[Rcpp::export]]
List radmix_fit_cpp(List data, int k, Nullable<NumericMatrix> init_q,
                    Nullable<NumericMatrix> init_p, bool fix_q,
                    bool fix_p, int max_iter, double rtol, int n_em_iter,
                    int n_qn, int n_cores, bool verbose,
                    Nullable<NumericMatrix> sup_q) {
  // The arguments have been validated in R; these checks only guard against
  // misuse from within the package. Every R error must happen here, before
  // any parallel region.
  GenoView geno(data);
  const size_t I = geno.nrow(), J = geno.ncol(), K = k;
  if (k < 1 || n_qn < 1 || n_cores < 1 || max_iter < 0 || n_em_iter < 0)
    stop("invalid arguments to radmix_fit_cpp()");
  if (I == 0 || J == 0)
    stop("There must be at least one individual and one locus.");
  // A training view of cross-validation precomputes which genotypes it
  // hides, once for all the passes of the fit (radmix_geno_access.h).
  geno.build_cv_mask(n_cores);
  const size_t nq = K * I, np = K * J;
  FitLog log = {verbose, {}};

  // Supervised analysis: the reference individuals (fixed[i] = 1) keep the
  // ancestry proportions given in sup_q (stored in q_sup, in the K x I
  // layout) throughout the fit: they are set before EM, the Q-side passes
  // copy them instead of updating them (radmix_kernels.h, "Fixed
  // individuals"), and they are restored after each quasi-Newton step. P is
  // estimated from all the individuals. fx is NULL in an unsupervised
  // analysis, which then runs exactly as without this option.
  std::vector<unsigned char> fixed;
  std::vector<double> q_sup;
  size_t n_fixed = 0;
  if (sup_q.isNotNull()) {
    NumericMatrix sq(sup_q.get());
    if ((size_t) sq.nrow() != I || (size_t) sq.ncol() != K)
      stop("invalid arguments to radmix_fit_cpp()");
    fixed.assign(I, 0);
    q_sup.assign(nq, 0.0);
    for (size_t i = 0; i < I; i++) {
      if (ISNAN(sq(i, 0))) continue;
      fixed[i] = 1;
      n_fixed++;
      for (size_t kk = 0; kk < K; kk++) q_sup[K * i + kk] = sq(i, kk);
    }
  }
  const unsigned char* fx = n_fixed > 0 ? fixed.data() : nullptr;
  // set the columns of the reference individuals of a q (K x I) to their
  // fixed values
  auto set_fixed = [&](std::vector<double>& qq) {
    for (size_t i = 0; i < I; i++)
      if (fixed[i])
        std::copy(&q_sup[K * i], &q_sup[K * (i + 1)], &qq[K * i]);
  };
  if (fx != nullptr) {
    log.add("Supervised analysis: " + std::to_string(n_fixed) +
            " reference individuals (fixed ancestry), " +
            std::to_string(I - n_fixed) + " individuals to estimate");
  }

  // the current estimate x0 = (q, p), and the estimates after one and two
  // sweeps of block updates, x1 = (q1, p1) and x2 = (q2, p2)
  std::vector<double> q(nq), p(np), q1(nq), p1(np), q2(nq), p2(np);

  // Random start: uniform values, each column of q rescaled to sum to one.
  // Drawn from R's random number generator (all of q, then all of p), so
  // set.seed() makes fits reproducible.
  for (size_t t = 0; t < nq; t++) q[t] = R::runif(0.0, 1.0);
  for (size_t i = 0; i < I; i++) {
    double sum = 0.0;
    for (size_t kk = 0; kk < K; kk++) sum += q[K * i + kk];
    for (size_t kk = 0; kk < K; kk++) q[K * i + kk] /= sum;
  }
  for (size_t t = 0; t < np; t++) p[t] = R::runif(0.0, 1.0);
  // the reference individuals start (and stay) at their fixed values, so
  // that EM estimates P from them (the random values drawn for them above
  // are discarded, which keeps the draws of the others in place)
  if (fx != nullptr) set_fixed(q);

  // EM priming, unless both blocks are given: a few EM steps move a random
  // start into a sensible region of the parameter space. Their
  // log-likelihoods are not needed (EM runs a fixed number of steps), and
  // each costs a pass over the genotypes, so they are only computed for a
  // verbose fit, every EM_LOG_EVERY steps.
  const int EM_LOG_EVERY = 5;
  // the log-likelihood of the current estimate, if the last EM step has
  // just computed it
  double ll_after_em = R_NaN;
  if (init_q.isNull() || init_p.isNull()) {
    for (int it = 1; it <= n_em_iter; it++) {
      em(geno, q.data(), p.data(), k, q1.data(), p1.data(), n_cores, fx);
      std::swap(q, q1);
      std::swap(p, p1);
      ll_after_em = R_NaN;
      if (verbose && it % EM_LOG_EVERY == 0) {
        ll_after_em = loglik(geno, q.data(), p.data(), k, n_cores);
        log.add("EM iteration " + std::to_string(it) +
                ": loglikelihood = " + fmt(ll_after_em, 12));
      }
      checkUserInterrupt();
    }
    if (n_em_iter > 0) {
      log.add("EM priming: " + std::to_string(n_em_iter) + " iterations");
    }
  }
  // Starting values given by the user replace their block (after EM, so that
  // when only one block is given the other one still starts from a sensible
  // value). Element (i, k) of the user's matrix is row k of column i of q.
  if (init_q.isNotNull()) {
    NumericMatrix iq(init_q.get());
    for (size_t i = 0; i < I; i++)
      for (size_t kk = 0; kk < K; kk++) q[K * i + kk] = iq(i, kk);
    // (the reference individuals keep their fixed values whatever init_q
    // says about them)
    if (fx != nullptr) set_fixed(q);
    ll_after_em = R_NaN;
  }
  if (init_p.isNotNull()) {
    NumericMatrix ip(init_p.get());
    for (size_t j = 0; j < J; j++)
      for (size_t kk = 0; kk < K; kk++) p[K * j + kk] = ip(j, kk);
    ll_after_em = R_NaN;
  }

  // Evaluating a point and starting the next iteration from it share a
  // pass: the first block update of an iteration, q1 = update of q given
  // (q, p), forms qp for every genotype, and the log-likelihood of (q, p)
  // costs little more on top (update_q_loglik()). So whenever the loop
  // needs the log-likelihood of a point from which the next iteration may
  // start (the start, the quasi-Newton point, x2), it computes q1 from that
  // point at the same time, speculatively; q1_ready records that q1 holds
  // the update of the current x0 = (q, p), so that the next iteration skips
  // its first pass. An accepted iteration then takes 4 passes over the
  // genotypes instead of 5, and a rejected one 5 instead of 6, with the same
  // estimates (q1 is identical either way) and the same log-likelihoods up
  // to rounding (the same terms, added up in another order). The update is
  // wasted only when no iteration follows (convergence, or a rejected
  // quasi-Newton point, whose update is replaced by that of x2), where the
  // pass costs a little more than the log-likelihood alone; it is not
  // computed when the iteration limit is reached, nor when q is fixed (there
  // is no update of q to compute), nor in a supervised analysis (fx): its
  // Q-side passes skip the reference individuals, so they cannot give the
  // log-likelihood, which a separate pass then computes.
  bool q1_ready = false;
  double ll_new;
  if (!ISNAN(ll_after_em)) {
    // (the last EM step of a verbose fit has already computed it)
    ll_new = ll_after_em;
  } else if (!fix_q && fx == nullptr && max_iter > 0) {
    ll_new = update_q_loglik(geno, q.data(), p.data(), k, q1.data(),
                             n_cores);
    q1_ready = true;
  } else {
    ll_new = loglik(geno, q.data(), p.data(), k, n_cores);
  }
  log.add("Initial loglikelihood: " + fmt(ll_new, 12));
  // the log-likelihood of the start, then at the end of every iteration
  std::vector<double> ll_history(1, ll_new);

  SecantPairs pairs(nq + np, n_qn);
  bool converged = false;
  int iter = 0;
  for (iter = 1; iter <= max_iter; iter++) {
    const double ll_prev = ll_new;

    // Two sweeps of block updates, x1 = F(x0) and x2 = F(x1). A sweep updates
    // all of q given p, then all of p given the new q, each column by a
    // constrained Newton step (radmix_kernels.cpp). A fixed block is copied,
    // and so are the columns of the reference individuals (fx). The first
    // update is usually already done (q1_ready, see above).
    if (fix_q) q1 = q;
    else if (!q1_ready)
      update_q(geno, q.data(), p.data(), k, q1.data(), n_cores, fx);
    q1_ready = false;
    if (fix_p) p1 = p;
    else update_p(geno, q1.data(), p.data(), k, p1.data(), n_cores);
    if (fix_q) q2 = q1;
    else update_q(geno, q1.data(), p1.data(), k, q2.data(), n_cores, fx);
    if (fix_p) p2 = p1;
    else update_p(geno, q2.data(), p1.data(), k, p2.data(), n_cores);

    // Quasi-Newton step (radmix_quasi_newton.cpp) from the last n_qn secant
    // pairs, kept in a circular buffer (fewer in the first iterations; their
    // order does not matter). The point overwrites x0, which is no longer
    // needed: if the step fails (a singular system), x0 is left untouched,
    // ll_qn stays at -Inf, and x2 is used instead.
    pairs.store((iter - 1) % n_qn, q.data(), p.data(), q1.data(), p1.data(),
                q2.data(), p2.data(), nq, np);
    const size_t n_used = std::min((size_t) iter, (size_t) n_qn);
    double ll_qn = R_NegInf;
    // whether to compute the next first update of q along with the
    // log-likelihood of the point that the iteration ends at (see above):
    // not after the last iteration allowed, nor with q fixed, nor with
    // reference individuals
    const bool speculate = !fix_q && fx == nullptr && iter < max_iter;
    if (qn_step(pairs.U.data(), pairs.V.data(), nq + np, n_used, q.data(),
                p.data(), q1.data(), p1.data(), k, I, J, q.data(), p.data(),
                n_cores)) {
      // A fixed block has u = v = 0, so the point keeps its value, but its
      // projection could move it by rounding errors: restore it exactly.
      if (fix_q) q = q1;
      if (fix_p) p = p1;
      // the same for the columns of the reference individuals
      if (fx != nullptr) set_fixed(q);
      // (q1 is no longer needed: x1 is only used by the secant pairs and the
      // quasi-Newton step, so it can receive the update of the point)
      if (speculate) {
        ll_qn = update_q_loglik(geno, q.data(), p.data(), k, q1.data(),
                                n_cores);
      } else {
        ll_qn = loglik(geno, q.data(), p.data(), k, n_cores);
      }
    }

    // The quasi-Newton point is kept if it improves on the log-likelihood of
    // x0; otherwise the iteration ends at x2. The comparison is with x0
    // rather than x2: a point slightly worse than x2 still makes progress,
    // and tends to pay off in the following iterations. It also means that
    // L(x2) is only needed when the point is rejected, saving a pass over
    // the genotypes in most iterations. (isfinite() rejects -Inf and NaN.)
    // If the point is rejected, the update of q computed with its
    // log-likelihood is discarded, and that of x2 is computed instead.
    double ll_basic = R_NaN;
    if (std::isfinite(ll_qn) && ll_prev < ll_qn) {
      ll_new = ll_qn;
      q1_ready = speculate;
    } else {
      if (speculate) {
        ll_basic = update_q_loglik(geno, q2.data(), p2.data(), k, q1.data(),
                                   n_cores);
      } else {
        ll_basic = loglik(geno, q2.data(), p2.data(), k, n_cores);
      }
      q1_ready = speculate;
      std::swap(q, q2);
      std::swap(p, p2);
      ll_new = ll_basic;
    }
    ll_history.push_back(ll_new);

    // convergence: the relative change in log-likelihood falls below rtol
    // (NaN, e.g. with no observed genotype, never converges)
    const double reldiff = std::fabs((ll_new - ll_prev) / ll_prev);
    log.add("Iteration " + std::to_string(iter) + ": loglikelihood = " +
            fmt(ll_new, 12) + ", reldiff = " + fmt(reldiff, 4) +
            " (basic = " +
            (ISNAN(ll_basic) ? std::string("not computed") :
             fmt(ll_basic, 12)) +
            ", QN = " + fmt(ll_qn, 12) + ")");
    if (reldiff < rtol) {
      converged = true;
      break;
    }
    checkUserInterrupt();
  }
  // (a loop that runs to the end leaves iter at max_iter + 1)
  if (iter > max_iter) iter = max_iter;
  if (converged) {
    log.add("Converged in " + std::to_string(iter) + " iterations.");
  } else {
    log.add("Failed to converge after " + std::to_string(iter) +
            " iterations.");
  }
  log.add("Loglikelihood: " + fmt(ll_new, 12));

  // the estimates in the user layout (individuals x K and loci x K)
  NumericMatrix q_out(I, K), p_out(J, K);
  for (size_t i = 0; i < I; i++)
    for (size_t kk = 0; kk < K; kk++) q_out(i, kk) = q[K * i + kk];
  for (size_t j = 0; j < J; j++)
    for (size_t kk = 0; kk < K; kk++) p_out(j, kk) = p[K * j + kk];
  return List::create(
    Named("q") = q_out, Named("p") = p_out, Named("loglik") = ll_new,
    Named("loglik_history") = wrap(ll_history), Named("n_iter") = iter,
    Named("converged") = converged, Named("log") = wrap(log.lines));
}
