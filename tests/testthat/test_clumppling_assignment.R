# Tests of the cluster alignment (src/clumppling_assignment.cpp)

# all permutations of 1:n
all_permutations <- function(n) {
  if (n == 1) {
    return(matrix(1L, 1, 1))
  }
  sub <- all_permutations(n - 1)
  do.call(rbind, lapply(seq_len(n), function(i) {
    cbind(i, matrix(setdiff(seq_len(n), i)[sub], nrow = nrow(sub)))
  }))
}

# all mappings of kq columns onto kp columns that cover every column
all_surjections <- function(kq, kp) {
  maps <- as.matrix(expand.grid(rep(list(seq_len(kp)), kq)))
  maps[apply(maps, 1, function(x) all(seq_len(kp) %in% x)), , drop = FALSE]
}

# the clumppling objective: sum of squared distances between matched columns,
# divided by N
alignment_obj <- function(p, q, idx) {
  sum(vapply(seq_len(ncol(q)), function(j) {
    sum((q[, j] - p[, idx[j]])^2)
  }, numeric(1))) / nrow(p)
}

test_that("the linear assignment solver finds the optimum", {
  withr::local_seed(1)
  for (n in 1:6) {
    perms <- all_permutations(n)
    for (rep in 1:5) {
      cost <- matrix(stats::runif(n * n), n, n)
      sol <- clumppling_lsap_cpp(cost)
      expect_setequal(sol, seq_len(n))
      brute <- min(apply(perms, 1, function(p) sum(cost[cbind(seq_len(n), p)])))
      expect_equal(sum(cost[cbind(seq_len(n), sol)]), brute)
    }
  }
  # integer costs with many ties
  cost <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2, 1), 3)
  expect_equal(sum(cost[cbind(1:3, clumppling_lsap_cpp(cost))]), 3)
  expect_error(clumppling_lsap_cpp(matrix(1, 2, 3)), "square")
  expect_error(clumppling_lsap_cpp(matrix(NaN, 2, 2)), "could not be solved")
})

test_that("align_ilp finds the optimal mapping, also with K_Q > K_P", {
  withr::local_seed(2)
  k_pairs <- list(c(2, 2), c(3, 3), c(4, 4), c(2, 3), c(3, 4), c(2, 4), c(3, 5))
  for (k in k_pairs) {
    for (rep in 1:3) {
      p <- clumppling_test_random_q(20, k[1])
      q <- clumppling_test_random_q(20, k[2])
      res <- clumppling_align_ilp_cpp(p, q)
      # every column of P receives at least one column of Q
      expect_setequal(unique(res$idx), seq_len(k[1]))
      expect_equal(res$obj, alignment_obj(p, q, res$idx))
      maps <- all_surjections(k[2], k[1])
      brute <- min(apply(maps, 1, function(idx) alignment_obj(p, q, idx)))
      expect_equal(res$obj, brute)
    }
  }
  # identical matrices are aligned to themselves
  p <- clumppling_test_random_q(30, 4)
  res <- clumppling_align_ilp_cpp(p, p[, c(3, 1, 4, 2)])
  expect_identical(res$idx, c(3L, 1L, 4L, 2L))
  expect_equal(res$obj, 0)
  expect_error(clumppling_align_ilp_cpp(p, p[-1, ]), "same number of rows")
  expect_error(clumppling_align_ilp_cpp(p, p[, 1:3]), "no fewer columns")
})

test_that("align_ilp_diff1 finds the best pair of clusters to merge", {
  withr::local_seed(3)
  for (kp in 2:4) {
    p <- clumppling_test_random_q(25, kp)
    q <- clumppling_test_random_q(25, kp + 1)
    res <- clumppling_align_ilp_diff1_cpp(p, q)
    # exactly one pair of columns of Q is merged
    expect_identical(sort(tabulate(res$idx, kp), decreasing = TRUE)[1:2], 2:1)
    # brute force: merge each pair, and align the merged matrix
    pairs <- utils::combn(kp + 1, 2)
    objs <- apply(pairs, 2, function(ab) {
      merged <- cbind(q[, -ab, drop = FALSE], q[, ab[1]] + q[, ab[2]])
      clumppling_align_ilp_cpp(p, merged)$obj
    })
    expect_equal(res$obj, min(objs))
  }
  expect_error(
    clumppling_align_ilp_diff1_cpp(p, cbind(q, q[, 1])),
    "exactly one more column"
  )
})

test_that("alignments match those of the python implementation", {
  ref <- clumppling_test_py_reference()$ilp
  inputs <- clumppling_test_unit_inputs()$ilp
  expect_length(ref, length(inputs))
  for (i in seq_along(inputs)) {
    p <- inputs[[i]]$P
    q <- inputs[[i]]$Q
    res <- clumppling_align_ilp_cpp(p, q)
    expect_identical(res$idx, as.integer(ref[[i]]$idx))
    expect_equal(res$obj, ref[[i]]$obj, tolerance = 1e-12)
    expect_equal(
      clumppling_align_q_wrt_p(p, q, res$idx, merge = TRUE)$q,
      ref[[i]]$merge_q,
      tolerance = 1e-14
    )
    if (ncol(q) - ncol(p) == 1) {
      res1 <- clumppling_align_ilp_diff1_cpp(p, q)
      expect_identical(res1$idx, as.integer(ref[[i]]$idx_diff1))
      expect_equal(res1$obj, ref[[i]]$obj_diff1, tolerance = 1e-12)
    }
    if (ncol(q) > ncol(p)) {
      no_merge <- clumppling_align_q_wrt_p(p, q, res$idx, merge = FALSE)
      expect_identical(no_merge$pattern, as.integer(ref[[i]]$nomerge_pattern))
      expect_equal(no_merge$q, ref[[i]]$nomerge_q, tolerance = 1e-14)
    }
  }
})

test_that("alignment within K returns consistent costs and patterns", {
  withr::local_seed(4)
  p <- clumppling_test_random_q(50, 3, alpha = 0.5)
  q_list <- list(
    p,
    p[, c(2, 3, 1)],
    clumppling_test_random_q(50, 3),
    p[, c(3, 2, 1)]
  )
  res <- clumppling_align_within_k(q_list)
  expect_identical(dim(res$cost), c(4L, 4L))
  expect_identical(dim(res$pattern), c(3L, 4L, 4L))
  expect_equal(res$cost, t(res$cost))
  expect_equal(diag(res$cost), rep(0, 4))
  expect_equal(res$cost[1, 2], 0)
  expect_equal(res$cost[2, 4], 0)
  expect_identical(res$pattern[, 1, 2], c(2L, 3L, 1L))
  for (i in 1:4) {
    expect_identical(res$pattern[, i, i], 1:3)
    for (j in setdiff(1:4, i)) {
      # the reverse alignment is the inverse permutation
      expect_identical(res$pattern[res$pattern[, i, j], j, i], 1:3)
      # the cost is the cost of the aligned matrices
      aligned <- clumppling_align_q_wrt_p(
        q_list[[i]], q_list[[j]], res$pattern[, i, j]
      )$q
      expect_equal(
        res$cost[i, j],
        clumppling_cost_membership(q_list[[i]], aligned)
      )
    }
  }
  # a single replicate
  single <- clumppling_align_within_k(q_list[1])
  expect_identical(single$cost, matrix(0, 1, 1))
  expect_identical(single$pattern, array(1:3, c(3, 1, 1)))
  expect_error(
    clumppling_align_within_k(list(p, p[-1, ])),
    "same dimensions"
  )
})
