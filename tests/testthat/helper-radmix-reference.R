# Pure R reference implementation of the kernels, on dense matrices.
# `g` is the decoded genotype matrix (individuals x loci, NA = missing),
# `n` the number of trials per individual, `q` is K x I and `p` is K x J.
# For genotype likelihoods, `g` is instead a list with the likelihoods e1
# and e2 (relative to the first genotype) and obs (FALSE for no data), all
# individuals x loci (see radmix_ref_gl_inputs()).

# the prepared dataset (`dat`, also available as `geno`, the name used for
# the argument of the kernels), and the decoded genotypes and numbers of
# trials used by the R references
radmix_ref_inputs <- function(gt) {
  dat <- radmix_prepare(gt)
  list(
    dat = dat,
    geno = dat,
    g = radmix_decode(dat),
    n = dat$n_trials
  )
}

# the same for a dataset of genotype likelihoods (see radmix_test_gl()),
# which keeps its likelihoods as matrices
radmix_ref_gl_inputs <- function(dat) {
  list(
    dat = dat,
    geno = dat,
    g = list(e1 = dat$e1, e2 = dat$e2, obs = dat$obs),
    n = rep(2, dat$n_indiv)
  )
}

# random starting values, in the layout of the kernels: q (K x I) with
# columns on the simplex and p (K x J) uniform in p_range; with seed = NULL,
# the global random number generator is used in its current state
radmix_ref_start <- function(n_indiv, n_loci, k, seed = 42,
                             p_range = c(0.05, 0.95)) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * n_loci, p_range[1], p_range[2]), k, n_loci)
  list(q = q, p = p)
}

# per-entry quantities shared by all the reference kernels
radmix_ref_terms <- function(g, n, q, p) {
  qp <- t(q) %*% p
  if (is.list(g)) {
    # genotype likelihoods: the expected counts of the alternate and
    # reference alleles (allele2 and allele1) given the likelihoods and qp
    w0 <- (1 - qp)^2
    w1 <- 2 * g$e1 * qp * (1 - qp)
    w2 <- g$e2 * qp^2
    f <- w0 + w1 + w2
    g0 <- ifelse(g$obs, (w1 + 2 * w2) / f, 0)
    h0 <- ifelse(g$obs, (2 * w0 + w1) / f, 0)
  } else {
    obs <- !is.na(g)
    g0 <- ifelse(obs, g, 0)
    h0 <- ifelse(obs, n - g0, 0)
  }
  # a = g / qp, b = (n - g) / (1 - qp), 0 when the coefficient is 0
  a <- ifelse(g0 > 0, g0 / qp, 0)
  b <- ifelse(h0 > 0, h0 / (1 - qp), 0)
  list(qp = qp, g0 = g0, h0 = h0, a = a, b = b)
}

radmix_ref_loglik <- function(g, n, q, p) {
  if (is.list(g)) {
    # genotype likelihoods, normalised to sum to 1
    qp <- t(q) %*% p
    f <- (1 - qp)^2 + 2 * g$e1 * qp * (1 - qp) + g$e2 * qp^2
    return(sum(ifelse(g$obs, log(f / (1 + g$e1 + g$e2)), 0)))
  }
  tt <- radmix_ref_terms(g, n, q, p)
  sum(ifelse(tt$g0 > 0, tt$g0 * log(tt$qp), 0)) +
    sum(ifelse(tt$h0 > 0, tt$h0 * log(1 - tt$qp), 0))
}

radmix_ref_em <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  # sum_j a p_kj + b (1 - p_kj), K x I
  sq <- p %*% t(tt$a) + (1 - p) %*% t(tt$b)
  # normalising by the column sums is the same as dividing by n_i times the
  # number of observed loci of individual i
  q_new <- q * sq
  q_new <- sweep(q_new, 2, colSums(q_new), "/")
  # expected counts of alternate (big_a) and reference (big_b) alleles from
  # each population at each locus (K x J); p is their ratio
  big_a <- p * (q %*% tt$a)
  big_b <- (1 - p) * (q %*% tt$b)
  p_new <- big_a / (big_a + big_b)
  list(q = q_new, p = p_new)
}

# gradient (K x I) of the log-likelihood with respect to each column of q,
# and the negative Hessian of each column (K x K x I), as returned by the
# kernels: with d1 = g / qp - (n - g) / (1 - qp) and
# d2 = g / qp^2 + (n - g) / (1 - qp)^2, grad_ki = sum_j d1_ij p_kj and
# hess_kk2i = sum_j d2_ij p_kj p_k2j
radmix_ref_grad_q <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  k <- nrow(q)
  grad <- p %*% t(tt$a - tt$b)
  w <- tt$a / tt$qp + tt$b / (1 - tt$qp)
  hess <- array(0, c(k, k, ncol(q)))
  for (i in seq_len(ncol(q))) {
    hess[, , i] <- p %*% (w[i, ] * t(p))
  }
  list(grad = grad, hess = hess)
}

# the same for the columns of p: gradient (K x J) and negative Hessian of
# each column (K x K x J)
radmix_ref_grad_p <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  k <- nrow(q)
  grad <- q %*% (tt$a - tt$b)
  w <- tt$a / tt$qp + tt$b / (1 - tt$qp)
  hess <- array(0, c(k, k, ncol(p)))
  for (j in seq_len(ncol(p))) {
    hess[, , j] <- q %*% (w[, j] * t(q))
  }
  list(grad = grad, hess = hess)
}

# compare every kernel with the pure R reference
radmix_expect_ref_match <- function(inp, k, tolerance = 1e-10) {
  st <- radmix_ref_start(inp$dat$n_indiv, inp$dat$n_loci, k)
  q <- st$q
  p <- st$p
  testthat::expect_equal(
    radmix_loglik_cpp(inp$geno, q, p, 1),
    radmix_ref_loglik(inp$g, inp$n, q, p),
    tolerance = tolerance
  )
  testthat::expect_equal(
    radmix_em_cpp(inp$geno, q, p, 1),
    radmix_ref_em(inp$g, inp$n, q, p),
    tolerance = tolerance
  )
  testthat::expect_equal(
    radmix_grad_q_cpp(inp$geno, q, p, 1),
    radmix_ref_grad_q(inp$g, inp$n, q, p),
    tolerance = tolerance
  )
  testthat::expect_equal(
    radmix_grad_p_cpp(inp$geno, q, p, 1),
    radmix_ref_grad_p(inp$g, inp$n, q, p),
    tolerance = tolerance
  )
}

# ALS on the factors of an SVD, as radmix_als_cpp() (src/radmix_als.cpp):
# za is loci x R and ua individuals x R, w the weights of the loci in the Q
# step, and init_p the starting allele frequencies (loci x K). Works in the
# user layout, solving one row at a time, but returns q (K x I) and p
# (K x J), as the C++ function.
radmix_ref_als <- function(za, ua, w, init_p, max_iter = 1000, tol = 1e-4) {
  clip <- function(a) pmin(pmax(a, 1e-5), 1 - 1e-5)
  map_q <- function(q) {
    q <- clip(q)
    q / rowSums(q)
  }
  # each row x_i of x minimising 0.5 x_i' a x_i - b_i' x_i in [0, 1]: the
  # unconstrained solution if it lies within the bounds (and a is well
  # conditioned), otherwise by cyclic coordinate descent from its current
  # value
  box_cd <- function(a, b, x) {
    direct <- rcond(a) > 1e-8
    if (direct) x_direct <- b %*% solve(a)
    for (i in seq_len(nrow(x))) {
      if (direct && all(x_direct[i, ] >= 0 & x_direct[i, ] <= 1)) {
        x[i, ] <- x_direct[i, ]
        next
      }
      for (s in 1:20) {
        change <- 0
        for (k in seq_len(ncol(x))) {
          new <- x[i, k] + (b[i, k] - sum(a[, k] * x[i, ])) / a[k, k]
          new <- min(max(new, 0), 1)
          change <- max(change, abs(new - x[i, k]))
          x[i, k] <- new
        }
        if (change < 1e-7) break
      }
    }
    x
  }
  q_step <- function(p, q) {
    wp <- p * w
    map_q(box_cd(crossprod(p, wp), ua %*% crossprod(za, wp), q))
  }
  p_step <- function(q, p) {
    clip(box_cd(crossprod(q), za %*% crossprod(ua, q), p))
  }
  p <- clip(init_p)
  k <- ncol(p)
  q <- q_step(p, matrix(1 / k, nrow(ua), k))
  q_prev <- q
  for (iter in seq_len(max_iter)) {
    q_use <- if (iter == 1) q else map_q(q + 0.8 * (q - q_prev))
    p <- p_step(q_use, p)
    q_new <- q_step(p, q)
    done <- sqrt(mean((q_new - q)^2)) < tol
    q_prev <- q
    q <- q_new
    if (done) break
  }
  # projected onto the feasible set of the fit, as radmix_als_cpp() does
  list(q = radmix_project_cpp(t(q), TRUE), p = t(p), n_iter = iter)
}
