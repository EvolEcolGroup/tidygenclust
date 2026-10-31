# Pure R reference implementation of the kernels, on dense matrices.
# `g` is the decoded genotype matrix (individuals x loci, NA = missing),
# `n` the number of trials per individual, `q` is K x I and `p` is K x J.

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

# random starting values with q on the simplex and p in (0, 1)
radmix_ref_start <- function(n_indiv, n_loci, k, seed = 42) {
  set.seed(seed)
  q <- matrix(stats::runif(k * n_indiv), k, n_indiv)
  q <- sweep(q, 2, colSums(q), "/")
  p <- matrix(stats::runif(k * n_loci, 0.05, 0.95), k, n_loci)
  list(q = q, p = p)
}

# per-entry quantities shared by all the reference kernels
radmix_ref_terms <- function(g, n, q, p) {
  qp <- t(q) %*% p
  obs <- !is.na(g)
  g0 <- ifelse(obs, g, 0)
  h0 <- ifelse(obs, n - g0, 0)
  # a = g / qp, b = (n - g) / (1 - qp), 0 when the coefficient is 0
  a <- ifelse(g0 > 0, g0 / qp, 0)
  b <- ifelse(h0 > 0, h0 / (1 - qp), 0)
  list(qp = qp, g0 = g0, h0 = h0, a = a, b = b)
}

radmix_ref_loglik <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  sum(ifelse(tt$g0 > 0, tt$g0 * log(tt$qp), 0)) +
    sum(ifelse(tt$h0 > 0, tt$h0 * log(1 - tt$qp), 0))
}

# the same, with 1 - qp written as sum_k q_k (1 - p_k): identical on the
# simplex, and the form whose derivatives OpenADMIXTURE.jl uses for the Q block
radmix_ref_loglik_homog <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  qomp <- t(q) %*% (1 - p)
  sum(ifelse(tt$g0 > 0, tt$g0 * log(tt$qp), 0)) +
    sum(ifelse(tt$h0 > 0, tt$h0 * log(qomp), 0))
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

# gradient (K x I) of the log-likelihood in its homogeneous form (see above)
# with respect to each column of q, and the negative Hessian of each column
# (K x K x I), as returned by the kernels
radmix_ref_grad_q <- function(g, n, q, p) {
  tt <- radmix_ref_terms(g, n, q, p)
  k <- nrow(q)
  grad <- p %*% t(tt$a) + (1 - p) %*% t(tt$b)
  wa <- tt$a / tt$qp
  wb <- tt$b / (1 - tt$qp)
  hess <- array(0, c(k, k, ncol(q)))
  for (i in seq_len(ncol(q))) {
    hess[, , i] <- p %*% (wa[i, ] * t(p)) +
      (1 - p) %*% (wb[i, ] * t(1 - p))
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
  st <- radmix_ref_start(nrow(inp$g), ncol(inp$g), k)
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
