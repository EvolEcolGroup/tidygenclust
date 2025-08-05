#' Install tools for `tidygenclust`
#'
#' `tidygenclust` relies on the `admxiture` and on python packages `fastmixture`
#' and `clumppling` for a number of
#' functionalities. We use `reticulate` to install them in conda
#' environments. As their dependencies are incompatible, we use two separate
#' conda environment, `ctidygenclust` (for `fastmixture` and `admixture`) and
#' `cclumppling` (for `clumppling`).
#' For each tool, default to the latest tested version of
#' these packages that have been tested to work with `tidegenclust`. It is
#' possible to provide a more recent github commit for a specific tool, but
#' this might lead to incompatibilities and errors.
#' @param reset a boolean used to reset the virtual environment. Only set
#' to TRUE if you have a broken virtual environment that you want to reset.
#' @param fastmixture_hash a string with the commit hash of the `fastmixture`
#' version to install. Default is the latest tested version.
#' @param clumppling_hash a string with the commit hash of the `clumppling`
#' version to install. Default is the latest tested version.
#' @returns NULL
#' @export

tgc_tools_install <-
  function(reset = FALSE,
           fastmixture_hash = "f913014669f4a235a1150669d4fbf0715bef42be",
           clumppling_hash = "a4bf351037fb569e2c2cb83c603a1931606d4d40") {
    if (.Platform$OS.type == "windows") {
      stop(
        "tidygenclust does not work on windows; use the Windows subsystem",
        "for Linux (WSL) instead"
      )
    }

    # check ctidygenclust does not exist
    if (reticulate::condaenv_exists("ctidygenclust")) {
      if (reset) {
        reticulate::conda_remove("ctidygenclust")
      } else {
        message(
          "The conda environment 'ctidygenclust' already exists. Use ",
          "'reset = TRUE' to reset it"
        )
        return(NULL)
      }
    }
    # check cclumppling does not exist
    if (reticulate::condaenv_exists("cclumppling")) {
      if (reset) {
        reticulate::conda_remove("cclumppling")
      } else {
        message(
          paste0(
            "The conda environment 'cclumppling' already exists. ",
            "Use 'reset = TRUE' to reset it"
          )
        )
        return(NULL)
      }
    }

    # install fastmixture
    reticulate::conda_create(
      envname = "ctidygenclust",
      packages = c("python>=3.10", "numpy>2.0.0", "cython>3.0.0"),
      channel = c("bioconda", "conda-forge", "defaults")
    )
    # create command line to install fastmixture
    fast_install_cmd <- paste0(
      "pip3 install ",
      "--upgrade --force-reinstall ",
      "git+https://github.com/Rosemeis/fastmixture.git@",
      fastmixture_hash
    )

    # on OSX we need to also install a suitable compiler for openmp
    if (Sys.info()["sysname"] == "Darwin") {
      reticulate::conda_install(
        envname = "ctidygenclust",
        packages = c("clang", "clangxx", "llvm-openmp"),
        channel = c("conda-forge", "defaults")
      )
      fast_install_cmd <- c(
        "export CC=clang",
        "export CXX=clang++",
        "export CFLAGS=-fopenmp",
        "export CXXFLAGS=-fopenmp",
        "export LDFLAGS=-fopenmp",
        fast_install_cmd
      )
    }

    ## https://github.com/rstudio/reticulate/issues/905
    reticulate::conda_run2(
      cmd_line = fast_install_cmd,
      envname = "ctidygenclust"
    )

    # if on osx or linux, install admixture
    if (.Platform$OS.type == "unix") {
      if ((Sys.info()["sysname"] == "Linux")) {
        # on linux we install admixture in ctidygenclust
        reticulate::conda_install(
          envname = "ctidygenclust",
          packages = c("admixture"),
          channel = c("bioconda")
        )
      } else if (Sys.info()["sysname"] == "Darwin") {
        # ADMIXTURE is only available for osx as x86 in bioconda
        # so we have to create a new environment and set it to x86_64
        reticulate::conda_create("cadmixture86",
          channel = c("bioconda", "conda-forge", "defaults")
        )
        reticulate::conda_run2(
          cmd = "conda",
          arg = "config --env --set subdir osx-64",
          envname = "cadmixture86"
        )
        # install admixture in the new environment
        reticulate::conda_install(
          envname = "cadmixture86",
          packages = c("admixture"),
          channel = c("bioconda")
        )
      }
    }



    #########################################################################
    # now install clumpling in its own conda environment
    # since its dependencies are not compatible with the ones of fastmixture
    reticulate::conda_create(
      envname = "cclumppling",
      packages = c("python==3.9", "numpy==1.24.0"),
      channel = c("bioconda", "conda-forge", "defaults")
    )
    # Install clumppling
    reticulate::conda_run2(
      cmd = "pip3",
      args = paste0(
        "install ",
        "--upgrade --force-reinstall ",
        "git+https://github.com/PopGenClustering/Clumppling.git@",
        clumppling_hash
      ),
      envname = "cclumppling"
    )

    #########################################################################
    # activate ctidygenclust with the python functions
    reticulate::use_condaenv("ctidygenclust", required = FALSE)
  }
