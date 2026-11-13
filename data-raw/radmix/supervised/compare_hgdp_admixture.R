# Supervised gt_radmix() against ADMIXTURE --supervised on HGDP pruned,
# K = 5: log-likelihood, Q of the individuals to estimate, iterations, time.
#
# Usage, from the root of the package:
#   Rscript data-raw/radmix/supervised/compare_hgdp_admixture.R \
#     <scratch dir> <case> <fit.rds>
# <scratch dir>/admix holds ADMIXTURE's run for <case> ("regions" or
# "pops"): hgdp_<case>.pop (written by hgdp_supervised.R, with "-" for the
# individuals to estimate), and the outputs of
#   admixture --supervised -s 1 hgdp_<case>.bed 5 -j1
# (hgdp_<case>.5.Q and its log, hgdp_<case>.5.log); <fit.rds> is the fit
# saved by hgdp_supervised.R for the same case.
args <- commandArgs(trailingOnly = TRUE)
scratch <- args[1]
case <- args[2]
fit <- readRDS(args[3])

adm_dir <- file.path(scratch, "admix")
labels <- readLines(file.path(adm_dir, paste0("hgdp_", case, ".pop")))
is_ref <- labels != "-"
adm_log <- readLines(file.path(adm_dir, paste0("hgdp_", case, ".5.log")))
q_adm <- as.matrix(utils::read.table(
  file.path(adm_dir, paste0("hgdp_", case, ".5.Q"))
))
ll_adm <- as.numeric(sub("Loglikelihood: ", "",
  grep("^Loglikelihood:", adm_log, value = TRUE)))
conv <- grep("^Converged in", adm_log, value = TRUE)
q_rad <- unname(unclass(fit$Q[[1]]))
# both programs order the columns by first appearance of the populations
stopifnot(identical(fit$ref_pops, unique(labels[is_ref])))
cat("case:", case, "\n")
cat("ADMIXTURE:", conv, "; loglik", format(ll_adm, digits = 15), "\n")
cat("gt_radmix:", grep("^Converged in", fit$log[[1]], value = TRUE),
  "; loglik", format(fit$loglik, digits = 15), "\n")
cat("loglik difference (gt_radmix - ADMIXTURE):", fit$loglik - ll_adm, "\n")
d <- abs(q_rad[!is_ref, ] - q_adm[!is_ref, ])
cat("Q of the", sum(!is_ref), "individuals to estimate: mean |diff| =",
  signif(mean(d), 3), ", max |diff| =", signif(max(d), 3), "\n")
cat("Q of the references: max |diff| =",
  max(abs(q_rad[is_ref, ] - q_adm[is_ref, ])), "\n")
