# Checks that the R replica of the fitting loop (fit_replica() in
# algo_lib.R) gives bit-identical iterates to radmix_fit(), for a random start
# and for a cross-validation refit, so that the trajectories it records are
# those of the package. Usage:
#   Rscript check_replica.R <lib> <dataset> <k> [n_cores]
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
dataset <- if (length(args) >= 2) args[2] else "hgdp_quick"
k <- if (length(args) >= 3) as.integer(args[3]) else 3L
n_cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
library(tidygenclust, lib.loc = lib)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(),
  value = TRUE
)))
source(file.path(script_dir, "algo_lib.R"))
x <- load_gt(dataset)
data <- tg$radmix_prepare(x)

# the package: a random start with seed 1
set.seed(1)
fit <- tg$radmix_fit(data, k = k, n_cores = n_cores)
# the replica: same start, 5 EM steps, then the loop
set.seed(1)
st <- random_start(k, data$n_indiv, data$n_loci)
for (i in 1:5) st <- tg$radmix_em_cpp(data, st$q, st$p, n_cores)
rep_fit <- fit_replica(data, st$q, st$p, n_cores = n_cores)
cat("main fit: iterations", fit$n_iter, rep_fit$n_iter,
  "| identical history:", identical(fit$loglik_history, rep_fit$loglik_history),
  "| identical Q:", identical(fit$q, t(rep_fit$q)),
  "| passes", fit_passes(fit), 5 + rep_fit$passes, "\n")

# cross-validation fold 1, warm-started from the fit
cv_seed <- 12345
train <- tg$radmix_set_cv(data, list(mode = "train", fold = 1, n_folds = 5,
  seed = cv_seed))
f1 <- tg$radmix_fit(train, k = k, init_q = fit$q, init_p = fit$p,
  n_cores = n_cores)
r1 <- fit_replica(train, t(fit$q), t(fit$p), n_cores = n_cores)
cat("fold 1: iterations", f1$n_iter, r1$n_iter, "| identical history:",
  identical(f1$loglik_history, r1$loglik_history), "| passes",
  fit_passes(f1, em_ran = FALSE), r1$passes, "\n")
