#' Tidy a `gt_clumppling` object
#'
#' A `tidy` method to extract information from a [`gt_clumppling`] object, and
#' return it as a tibble. It can extract:
#' - 'modes': all the modes detected by [gt_clumppling()]. The modes have label
#' 'KxMy', where 'x' and 'y' represent the K value and the mode rank.
#' - 'major_modes': modes of rank 1 for each K.
#' - 'Q_modes', 'q_modes': a list of q matrices, one per mode, each tidied
#' into a tibble
#' - 'Q_major_modes', 'q_major_modes': the same output as 'Q_modes' but
#' subsetted to only the
#' major modes.
#'
#' @param x the [`gt_clumppling`] object
#' @param matrix a string defining the information to be extracted, one of:
#'   "modes", "major_modes", "Q_modes" (or "q_modes"), "Q_major_modes" (or
#'   "q_major_modes").
#' @param ... Additional arguments. Not used. Needed to match generic signature
#'   only.
#' @returns a [tibble::tibble] of the information of interest (for "Q_modes"
#'   and "Q_major_modes", a named list of tibbles, one per mode)
#' @rdname tidy_gt_clumppling
#' @export

tidy.gt_clumppling <- function(
  x,
  matrix = c(
    "modes",
    "major_modes",
    "Q_modes",
    "q_modes",
    "Q_major_modes",
    "q_major_modes"
  ),
  ...
) {
  rlang::check_dots_empty()
  matrix <- match.arg(matrix)
  if (matrix == "modes") {
    # the K and the rank of each mode are parsed from its label (e.g. "K3M2")
    all_modes <- names(x$aligned_modes)
    tibble::tibble(
      k = as.numeric(gsub(".*K(\\d+)M.*", "\\1", all_modes)),
      m = as.numeric(gsub(".*M(\\d+).*", "\\1", all_modes)),
      label = all_modes
    )
  } else if (matrix == "major_modes") {
    tidy(x, matrix = "modes") %>% dplyr::filter(.data$m == 1)
  } else if (matrix == "Q_modes" || matrix == "q_modes") {
    # TODO we need individual ids, either as row names or we add arbitrary ones
    lapply(x$aligned_modes, tidy_q)
  } else if (matrix == "Q_major_modes" || matrix == "q_major_modes") {
    lapply(x$aligned_modes[tidy(x, matrix = "major_modes")$label], tidy_q)
  }
}

# Tidy the Q matrix of a mode (a data.frame with columns V1, V2, ..., one per
# cluster) into a long tibble with one row per individual and cluster, and
# columns 'id' (the index of the individual, taken from the row names, so
# that it is kept when the object is subsetted by individuals), 'q' (the
# cluster number, as a string) and 'percentage' (the membership)
tidy_q <- function(x) {
  q_tbl <- x %>%
    tibble::as_tibble() %>%
    # rename V1, V2, ... to .Q1, .Q2, ...
    dplyr::rename_with(~ sub("^V", ".Q", .x)) %>%
    # add the pops data for plotting here if needed
    dplyr::mutate(id = as.integer(rownames(x))) %>%
    # TODO we could add here the optional group info
    tidyr::pivot_longer(
      cols = dplyr::starts_with(".Q"),
      names_to = "q",
      values_to = "percentage"
    ) %>%
    # keep only the cluster number (e.g. ".Q2" becomes "2")
    dplyr::mutate(q = sub("^\\.Q", "", .data$q))
  q_tbl
}
