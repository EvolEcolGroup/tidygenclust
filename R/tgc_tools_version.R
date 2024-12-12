#' Print version of the python tools installed by tidygenclust
#' @description This function prints the version of the python tools installed by tidygenclust
#' @return A character vector with the version of the python tools installed by tidygenclust
#' @export

tgc_tools_version <- function() {

  if (.Platform$OS.type %in% c("unix", "darwin")) {
    reticulate::conda_run2(cmd = "admixture", args = c("--version"), envname = "ctidygenclust",
                         echo=FALSE)
  }
  reticulate::conda_run2(cmd = "fastmixture", args = c("--version"), envname = "ctidygenclust",
                         echo=FALSE)
  cat("cclumpling ")
  reticulate::conda_run2(cmd = "pip", args = c("show clumppling | grep Version"), envname = "cclumppling",
                         echo=FALSE)
}
