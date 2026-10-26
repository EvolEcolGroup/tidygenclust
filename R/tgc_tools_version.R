#' Print version of the optional tools installed by `tidygenclust`
#' @description This function returns the version of the optional tools
#'   installed by [tgc_tools_install()]: `fastmixture` and `ADMIXTURE`. Tools
#'   that are not installed have version NA.
#' @return A list with the version of the optional tools installed by
#' `tidygenclust`
#' @export

tgc_tools_version <- function() {
  # run a command in the conda environment of a tool, returning NA if the
  # tool is not installed
  tool_output <- function(tool, cmd, args) {
    if (!tgc_tool_installed(tool)) {
      return(NA)
    }
    # the output lines of the command (a non-zero exit status gives a
    # warning, which is also turned into NA)
    tryCatch(
      reticulate::conda_run2(
        cmd = cmd,
        args = args,
        envname = tgc_tool_env(tool),
        echo = FALSE,
        intern = TRUE
      ),
      error = function(e) NA,
      warning = function(w) NA
    )
  }

  fast_ver <- tool_output("fastmixture", "fastmixture", "--version")
  # ADMIXTURE prints a banner before its version number, which is the last
  # line of the output
  adm_ver <- tool_output("admixture", "admixture", "--version")
  list(
    fastmixture = fast_ver,
    admixture = adm_ver[length(adm_ver)]
  )
}
