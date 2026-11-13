# Supervised gt_radmix() against ADMIXTURE --supervised on EUR_subset, K = 3
#
# EUR_subset (379 individuals, 54,051 SNPs) has no population labels, so the
# reference individuals are chosen with a documented rule: from the stored
# unsupervised ADMIXTURE fit at K = 3 (tests/testthat/fixtures/
# admixture_eur_reference.rds), every individual with q > 0.9 in a cluster
# is a reference of that cluster (15, 59 and 28 individuals); the other 277
# are estimated. The clusters are labelled "mid", "zeta" and "alpha"
# (clusters 1, 2 and 3), so that the order of first appearance in the file
# (zeta, alpha, mid) differs from the alphabetical one, which shows how
# ADMIXTURE orders its columns.
#
# Usage, from the root of the package:
#   Rscript data-raw/radmix/supervised/compare_eur.R <library> <scratch dir> \
#     [n_cores]
# The scratch dir receives the .pop file, the ADMIXTURE outputs and the
# gen_tibble backing file. ADMIXTURE is run only if its output is missing.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scratch <- args[2]
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 1L
library(tidygenclust, lib.loc = lib)

eur_bed <- "/home/andrea/.cache/R/tidygenclust/EUR_subset.bed"
admixture <- "/home/andrea/.local/bin/admixture"
fixture <- readRDS("tests/testthat/fixtures/admixture_eur_reference.rds")

# the reference labels (NA for the individuals to estimate)
q3 <- fixture$runs[[2]]$q_millionths / 1e6
cluster <- apply(q3, 1, which.max)
is_ref <- apply(q3, 1, max) > 0.9
labels <- c("mid", "zeta", "alpha")[cluster]
labels[!is_ref] <- NA

# ADMIXTURE --supervised (1 thread), with the .pop file next to the bed
adm_dir <- file.path(scratch, "admix_eur")
dir.create(adm_dir, showWarnings = FALSE, recursive = TRUE)
for (ext in c("bed", "bim", "fam")) {
  file.symlink(sub("bed$", ext, eur_bed),
    file.path(adm_dir, paste0("eur.", ext)))
}
writeLines(ifelse(is.na(labels), "-", labels), file.path(adm_dir, "eur.pop"))
if (!file.exists(file.path(adm_dir, "eur.3.Q"))) {
  old <- setwd(adm_dir)
  t_adm <- system.time(system2(admixture,
    c("--supervised", "-s", "1", "eur.bed", "3", paste0("-j", n_cores)),
    stdout = "eur.3.log", stderr = "eur.3.log"
  ))
  writeLines(format(t_adm[["elapsed"]]), "eur.3.time")
  setwd(old)
}
adm_log <- readLines(file.path(adm_dir, "eur.3.log"))
adm_q <- as.matrix(utils::read.table(file.path(adm_dir, "eur.3.Q")))
adm_ll <- as.numeric(sub("Loglikelihood: ", "",
  grep("^Loglikelihood:", adm_log, value = TRUE)))
adm_iter <- as.integer(sub("Converged in ([0-9]+) iterations.*", "\\1",
  grep("^Converged in", adm_log, value = TRUE)))
adm_time <- as.numeric(readLines(file.path(adm_dir, "eur.3.time")))

eur <- tidypopgen::gen_tibble(eur_bed,
  backingfile = tempfile("eur_", tmpdir = scratch), quiet = TRUE
)
t_rad <- system.time(fit <- gt_radmix(eur,
  k = 3, supervised = labels, seed = 1, n_cores = n_cores
))[["elapsed"]]
q_rad <- unclass(fit$Q[[1]])

cat("ADMIXTURE column order: ",
  paste(unique(labels[!is.na(labels)]), collapse = ", "),
  " (first appearance); gt_radmix ref_pops: ",
  paste(fit$ref_pops, collapse = ", "), "\n", sep = "")
res <- data.frame(
  program = c("ADMIXTURE", "gt_radmix"),
  loglik = c(adm_ll, fit$loglik),
  iterations = c(adm_iter, as.integer(sub(".* in ([0-9]+) iterations.*",
    "\\1", grep("^Converged in", fit$log[[1]], value = TRUE)))),
  time_s = c(adm_time, t_rad),
  n_cores = n_cores
)
print(res, digits = 12)
unk <- is.na(labels)
cat("Q of the individuals to estimate: max |diff| =",
  max(abs(q_rad[unk, ] - adm_q[unk, ])), ", mean |diff| =",
  mean(abs(q_rad[unk, ] - adm_q[unk, ])), "\n")
cat("Q of the references: max |diff| =",
  max(abs(q_rad[!unk, ] - adm_q[!unk, ])), "\n")
cat("Our log-likelihood of ADMIXTURE's (rounded) estimates is not",
  "computed here: ADMIXTURE's P is not read\n")
utils::write.csv(res, file.path(scratch, "compare_eur.csv"), row.names = FALSE)
