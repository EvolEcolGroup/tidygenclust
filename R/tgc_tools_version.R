#' Print version of the python tools installed by tidygenclust
#' @description This function prints the version of the python tools installed by tidygenclust
#' @return A character vector with the version of the python tools installed by tidygenclust
#' @export

tgc_tools_version <- function() {
  reticulate::conda_run2(cmd = "fastmixture", args = c("--version"), envname = "ctidygenclust",
                         echo=FALSE)
  cat("clumpling ")
  reticulate::conda_run2(cmd = "pip", args = c("show clumppling | grep Version"), envname = "ctidygenclust",
                         echo=FALSE)
}
