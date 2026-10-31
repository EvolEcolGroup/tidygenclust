# Create the example gt_admix objects in inst/extdata from the Q files of the
# 'capeverde' (ADMIXTURE, 200 replicates, K = 2 to 5) and 'chicken_gapK'
# (STRUCTURE, 60 replicates, K = 17, 19 and 21) examples of clumppling
# (data-raw/clumppling/capeverde.zip and chicken_gapk.zip, which previous
# versions of tidygenclust distributed in inst/extdata).
#
# Requirements: tidypopgen.
# Run from the root of the package:
#   Rscript data-raw/clumppling/create_example_gt_admix.R
#
# Files are added in alphabetical order (as python sorts them), which is the
# order in which the python reference results
# (tests/testthat/fixtures/clumppling_python_reference.rds) were computed, so
# that replicates are numbered in the same way.

# unzip one of the example archives, and return the sorted paths of its files
unzip_example <- function(file) {
  exdir <- tempfile()
  utils::unzip(file.path("data-raw", "clumppling", file), exdir = exdir)
  sort(list.files(exdir, full.names = TRUE), method = "radix")
}

# parse the lines of an .indivq file, or of the inferred ancestry of a
# STRUCTURE output file, which have the columns: index, label, % missing,
# population, ':' and memberships
parse_q_lines <- function(lines) {
  tab <- utils::read.table(
    text = lines,
    colClasses = c(V2 = "character", V4 = "character")
  )
  stopifnot(all(tab$V5 == ":"))
  list(
    label = tab$V2,
    pop = tab$V4,
    q = as.matrix(tab[, -(1:5)])
  )
}

# the lines with the inferred ancestry of individuals in a STRUCTURE output
# file (between the header of the section and the following empty line)
structure_lines <- function(file) {
  lines <- readLines(file, warn = FALSE)
  start <- match("Inferred ancestry of individuals:", lines) + 2
  n_lines <- match("", lines[-seq_len(start - 1)]) - 1
  lines[start - 1 + seq_len(n_lines)]
}

# a gt_admix object, with the same layout as those created by tidypopgen
as_gt_admix <- function(runs, id, group, algorithm) {
  q_list <- lapply(runs, function(x) tidypopgen::q_matrix(x$q))
  x <- list(
    k = vapply(q_list, ncol, integer(1)),
    Q = q_list,
    id = id,
    group = group,
    algorithm = algorithm
  )
  class(x) <- c("gt_admix", "list")
  x
}

# capeverde: the individual id is the label, and the population is numbered
# in the files (the names are in the capeverde_pops dataset)
load(file.path("data", "capeverde_pops.rda"))
capeverde_runs <- lapply(unzip_example("capeverde.zip"), function(file) {
  parse_q_lines(readLines(file, warn = FALSE))
})
pops <- capeverde_runs[[1]]$pop
stopifnot(
  all(vapply(capeverde_runs, function(x) identical(x$pop, pops), logical(1))),
  identical(
    rle(pops)$lengths,
    rle(capeverde_pops)$lengths
  )
)
capeverde <- as_gt_admix(capeverde_runs,
  id = capeverde_runs[[1]]$label,
  group = capeverde_pops,
  algorithm = "ADMIXTURE"
)
saveRDS(capeverde, "inst/extdata/capeverde_admix.rds", compress = "xz")

# chicken: STRUCTURE rounds memberships to 3 decimals, so rows do not sum to
# 1; they are rescaled, which gives the same results as clumppling on the
# original files (which it validates with a tolerance of 1e-2) to within
# 1e-10
chicken_runs <- lapply(unzip_example("chicken_gapk.zip"), function(file) {
  x <- parse_q_lines(structure_lines(file))
  x$q <- x$q / rowSums(x$q)
  x
})
labels <- chicken_runs[[1]]$label
stopifnot(all(vapply(chicken_runs, function(x) {
  identical(x$label, labels)
}, logical(1))))
chicken <- as_gt_admix(chicken_runs,
  id = labels,
  group = chicken_runs[[1]]$pop,
  algorithm = "STRUCTURE"
)
saveRDS(chicken, "inst/extdata/chicken_admix.rds", compress = "xz")
