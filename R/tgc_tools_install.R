#' Install tools for `tidygenclust`
#'
#' `tidygenclust` relies on the python packages `fastmixture` and `clumppling` for a number of
#' functionalities. We use `reticulate` to install them in a conda
#' environment. We also install the other tools in the same environment to
#' make sure that we have everything in the same place.  By default, we use a
#' conda environment called 'ctidygenclust'. Unless
#' you know what you are doing, we recommend that you use this default. For
#' more details, see \href{https://rstudio.github.io/reticulate/articles/python_dependencies.html}{here}.
#'
#' @param envname the name of the virtual environment to use (defaults to 'ctidygenclust')
#' @param reset a boolean used to reset the virtual environment. Only set
#' to TRUE if you have a broken virtual environment that you want to reset.
#' @returns NULL
#' @export

tgc_tools_install <- function(envname = "ctidygenclust",
                             reset = FALSE) {
  if (reticulate::condaenv_exists(envname)){
    if (reset){
      reticulate::conda_remove(envname)
    } else {
      message("The conda environment ", envname, " already exists. Use 'reset = TRUE' to reset it")
      return(NULL)
    }
  }

  reticulate::conda_create(envname = "ctidygenclust",
                           packages = c("python=3.11", "numpy", "cython"),
                           channels = c("defaults", "bioconda"))
  ## https://github.com/rstudio/reticulate/issues/905
  reticulate:::conda_run2(cmd = "pip3",
                          args = paste0("install ",
                                        "--upgrade --force-reinstall ",
                                        "git+https://github.com/Rosemeis/fastmixture.git@main"),
                          envname = envname)

  # Install clumppling
  reticulate:::conda_run2(cmd = "pip3",
                          args = paste0("install ",
                                        "--upgrade --force-reinstall ",
                                        "git+https://github.com/PopGenClustering/Clumppling.git@master"),
                          envname = envname)

  # activate it
  reticulate::use_condaenv(envname, required = FALSE)
}
