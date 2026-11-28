# Helpers shared by several functions of the package: checks of the
# arguments, with errors that name the argument, and the group of each
# individual of a gen_tibble for gt_admix objects

# Check that an argument is a whole number of at least `min` (or, with
# scalar = FALSE, a vector of them), with an error naming the argument
check_count <- function(x, name, min = 1, scalar = TRUE, message = NULL) {
  ok <- is.numeric(x) && length(x) >= 1 && (!scalar || length(x) == 1) &&
    all(is.finite(x)) && all(x == round(x)) && all(x >= min)
  if (!ok) {
    if (is.null(message)) {
      message <- paste0("'", name, "' should be an integer of at least ", min)
    }
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

# Check that an argument is a single TRUE or FALSE
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("'", name, "' should be TRUE or FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

# The group of each individual of a grouped gen_tibble, or of any other
# grouped tibble (e.g. genotype likelihoods from another package), NULL if
# it is not grouped, for the `group` element of a gt_admix object. As in
# gt_admixture(), a single grouping variable is supported
gt_single_group <- function(x) {
  if (!inherits(x, "grouped_df")) {
    return(NULL)
  }
  group_vars <- dplyr::group_vars(x)
  if (length(group_vars) != 1) {
    stop(
      "'x' should be grouped by a single variable (it is grouped by ",
      paste(group_vars, collapse = ", "), ")",
      call. = FALSE
    )
  }
  x[[group_vars]]
}
