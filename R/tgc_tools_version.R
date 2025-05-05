#' Print version of the python tools installed by tidygenclust
#' @description This function prints the version of the python tools installed by tidygenclust
#' @return A list with the version of the python tools installed by tidygenclust
#' @export

tgc_tools_version <- function() {
  fast_ver <- reticulate::conda_run2(cmd = "fastmixture", args = c("--version"), envname = "ctidygenclust",
                         echo=FALSE, intern = TRUE)
  clump_ver<- reticulate::conda_run2(cmd = "pip", args = c("show clumppling | grep Version"), envname = "cclumppling",
                         echo=FALSE, intern = TRUE)
  tgc_vers <- list(fastmixture = fast_ver, clumppling = clump_ver)
  if (.Platform$OS.type %in% c("unix", "darwin")) {
    adm_ver <- reticulate::conda_run2(cmd = "admixture", args = c("--version"), envname = "ctidygenclust",
                                      echo=FALSE, intern=TRUE)
    tgc_vers[["admxiture"]] <- adm_ver[length(adm_ver)]
  }
  return(tgc_vers)
}
