#' Install tools for `tidygenclust`
#'
#' `tidygenclust` relies on `ADMIXTURE` and on python packages `fastmixture` and
#' `clumppling` for a number of functionalities. We use `reticulate` to install
#' them in conda environments. As their dependencies are incompatible, we use
#' two separate conda environments, `ctidygenclust` (for `fastmixture` and
#' `admixture`) and `cclumppling` (for `clumppling`). Additionally, for silicon
#' Macs, `ADMIXTURE` is installed in a separate conda environment
#' `cadmixture86`, as it is only available for OSX as x86 in bioconda.
#' @details For each tool, default to the latest tested version of these
#'   packages that have been tested to work with `tidegenclust`. It is possible
#'   to provide a more recent github commit for a specific tool, but this might
#'   lead to incompatibilities and errors.
#'
#'   We have found installation on OSX to be tricky, so we provide two methods
#'   for installing `fastmixture` on OSX: `reticulate` and `conda_yaml`. The
#'   `reticulate` method uses the `reticulate::conda_run2()` function to run
#'   installation commands, while the `conda_yaml` method creates a conda
#'   environment directly with conda. If the `reticulate` method fails, you can
#'   use the `conda_yaml` method to create the environment directly with conda.
#'   For OSX, you might also need to install a suitable compiler for openmp
#'   using `brew` in `bash`, setting the correct paths to use it:
#'
#'   `brew install llvm libomp`
#'
#' @param reset a boolean used to reset the virtual environment. Only set to
#'   TRUE if you have a broken virtual environment that you want to reset.
#' @param fastmixture_hash a string with the commit hash of the `fastmixture`
#'   version to install. Default is the latest tested version.
#' @param clumppling_hash a string with the commit hash of the `clumppling`
#'   version to install. Default is the latest tested version.
#' @param conda_yaml A vector of string with the path for the `yaml` files used
#'   for installation (one for fastmixture, one for clumppling, and, for osx
#'   only, one for admixture in 86 compatibility model). The value of each
#'   string can either be "auto" (the default, select the latest available yaml
#'   that has been tested to work on a given operating system), NULL (let conda
#'   attempt to install the latest available packages from the internet), or a
#'   path to a specific conda yaml found in "inst/env_snapshots". For "auto" or
#'   NULL, it is also possible to provide a vector of lenght one, and the value
#'   will apply to all the conda environments (e.g. `conda_yaml="auto"` is
#'   equivalent to  `conda_yaml = c("auto", "auto")` in linux and
#'   `conda_yaml = c("auto", "auto", "auto")`.
#' @param ci_install a boolean indicating if the installation is being run on
#'   continuous integration (CI) services. Default is FALSE. If TRUE, the
#'   function will look for the conda yaml file in the `inst/python` folder of
#'   the package source directory, rather than in the installed package
#'   directory. This is useful when testing the package on CI services.
#' @returns NULL
#' @export

tgc_tools_install <-
  function(reset = FALSE,
           fastmixture_hash = "29e04339ce6ddf750ee4e06f8aabe40335e0d0ee",
           clumppling_hash = "2d24e0b2f6ddfcb51a436df96a06d5f57d18d20a",
           conda_yaml = "auto",
           ci_install = FALSE) {
    # give error for windows
    if (.Platform$OS.type == "windows") {
      stop(
        "tidygenclust does not work on windows; use the Windows subsystem",
        "for Linux (WSL) instead"
      )
    }

    # if conda_yaml is a vector of length 1, repeat it to match the number of
    # conda environments
    if (length(conda_yaml) == 1) {
      if (Sys.info()["sysname"] == "Darwin") {
        conda_yaml <- rep(conda_yaml, 3)
      } else {
        conda_yaml <- rep(conda_yaml, 2)
      }
    }

    # now, if any value of conda yaml is "auto", replace it with the path to the
    # latest available yaml
    for (i in seq(along = conda_yaml)) {
      if (conda_yaml[i] == "auto") {
        if (Sys.info()["sysname"] == "Darwin") {
          if (i == 1) {
            conda_yaml[i] <- tail(
              list.files(
                system.file("env_snapshots/", package = "tidygenclust"),
                pattern = "ctidygenclust_osx_",
                full.names = TRUE
              ), 1
            )
          } else if (i == 2) {
            conda_yaml[i] <- tail(
              list.files(
                system.file("env_snapshots/", package = "tidygenclust"),
                pattern = "cclumppling_osx_",
                full.names = TRUE
              ), 1
            )
          } else if (i == 3) {
            conda_yaml[i] <- tail(
              list.files(
                system.file("env_snapshots/", package = "tidygenclust"),
                pattern = "cadmixture86_osx_",
                full.names = TRUE
              ), 1
            )
          }  
        } else { # if we are on linux, we only have two conda environments
          if (i == 1) {
            conda_yaml[i] <- tail(
              list.files(
                system.file("env_snapshots/", package = "tidygenclust"),
                pattern = "ctidygenclust_linux_",
                full.names = TRUE
              ), 1
            )
          } else if (i == 2) {
            conda_yaml[i] <- tail(
              list.files(
                system.file("env_snapshots/", package = "tidygenclust"),
                pattern = "cclumppling_linux_",
                full.names = TRUE
              ), 1
            )
          }
        }
      }
    }

    # if the values in conda_yaml are not null, check that the files exist
    if (!is.null(conda_yaml)) {
      for (i in seq(along = conda_yaml)) {
        if (!file.exists(conda_yaml[i])) {
          stop(
            paste0(
              "The conda yaml file '",
              conda_yaml[i],
              "' does not exist. Please provide a valid path to a conda yaml file."
            )
          )
        }
      }
    }

    # output the names of the yaml files being used for installation
    message(
      paste0(
        "Using the following conda yaml files for installation: ",
        paste(conda_yaml, collapse = ", ")
      )
    )
    
    
    
    # for osx, check that we have installed the right packages in brew
    if (Sys.info()["sysname"] == "Darwin") {
      brew_pkgs <- c("llvm", "libomp")
      if (!all(brew_installed(brew_pkgs))) {
        stop(
          paste0(
            "On OSX, please install the following packages with brew: ",
            paste(brew_pkgs, collapse = ", "),
            ". See ?tgc_tools_install for details"
          )
        )
      }
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
    # check that cadmixture86 does not exist
    if (reticulate::condaenv_exists("cadmixture86")) {
      if (reset) {
        reticulate::conda_remove("cadmixture86")
      } else {
        message(
          "The conda environment 'cadmixture86' already exists. Use ",
          "'reset = TRUE' to reset it"
        )
        return(NULL)
      }
    }

    ############################################
    # install fastmixture
    ############################################
    
      # create a conda environment with the necessary packages
      reticulate::conda_create(
        envname = "ctidygenclust",
        packages = c("python==3.11", "numpy>2.0.0", "cython>3.0.0", "pip"),
        channel = c("bioconda", "conda-forge", "defaults"),
        environment = conda_yaml[1]
      )
      # create command line to install fastmixture
      fast_install_cmd <- paste0(
        "python -m pip install ",
        "--upgrade --force-reinstall ",
        "git+https://github.com/Rosemeis/fastmixture.git@",
        fastmixture_hash
      )

      # on OSX we need to also install a suitable compiler for openmp
      if (Sys.info()["sysname"] == "Darwin") {
        # install clang and llvm-openmp
        reticulate::conda_install(
          envname = "ctidygenclust",
          packages = c("clang", "clangxx", "llvm-openmp"),
          channel = c("conda-forge", "defaults")
        )
        fast_install_cmd <- c(
          "export PATH=\"/opt/homebrew/opt/llvm/bin:$PATH\"",
          "export CC=clang",
          "export CXX=clang++",
          "export CFLAGS=-fopenmp",
          "export CXXFLAGS=-fopenmp",
          "export LDFLAGS=-fopenmp",
          fast_install_cmd
        )
      }

      # now install fastmixture in the conda environment
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
          channel = c("bioconda", "conda-forge", "defaults"),
          environment = conda_yaml[3]
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
    # now install clumppling in its own conda environment
    # since its dependencies are not compatible with the ones of fastmixture

    # find out operating system and release
    release <- Sys.info()[["release"]]
    version_major <- as.numeric(strsplit(release, "\\.")[[1]][1])

    if (Sys.info()[["sysname"]] == "Darwin" && version_major <= 23) {
      # macOS 14 or older
      reticulate::conda_create(
        envname = "cclumppling",
        packages = c("python==3.9", "pip"),
        channel = c("bioconda", "conda-forge", "defaults")
      )
      # NOTE we don't use the yaml for this old version of osx; we rely on a few
      # hacks to install the right version of setuptools and cvxopt, which are
      # not compatible with the latest version of python.
      reticulate::conda_run2(
        cmd = "python",
        args = paste0(
          "-m pip ",
          "install ",
          "--upgrade --force-reinstall ",
          "setuptools==80.10.2 ",
          "cvxopt==1.3.2 ",
          "git+https://github.com/PopGenClustering/Clumppling.git@",
          clumppling_hash
        ),
        envname = "cclumppling"
      )
    } else {
      # all other operating systems and macOS 15 or newer
      reticulate::conda_create(
        envname = "cclumppling",
        packages = c("python==3.12", "pip"),
        channel = c("bioconda", "conda-forge", "defaults"),
        environment = conda_yaml[2]
      )
      # Install clumppling
      reticulate::conda_run2(
        cmd = "python",
        args = paste0(
          "-m pip ",
          "install ",
          "--upgrade --force-reinstall ",
          "git+https://github.com/PopGenClustering/Clumppling.git@",
          clumppling_hash
        ),
        envname = "cclumppling"
      )
    }

    # check clumppling has successfully installed and warn user if it hasn't
    out <- system2(
      command = "conda",
      args = c("run", "-n", "cclumppling", "conda", "list", "clumppling"),
      stdout = FALSE,
      stderr = FALSE
    )

    if (out != 0) {
      warning(paste0(
        "clumppling has not been succesfully installed ",
        "in your conda environment"
      ))
    }

    #########################################################################
    # activate ctidygenclust with the python functions
    reticulate::use_condaenv("ctidygenclust", required = TRUE)
  }

# check if package is installed with brew
brew_installed <- function(pkgs) {
  vapply(pkgs, function(pkg) {
    status <- system2("brew",
      args = c("list", "--versions", pkg),
      stdout = FALSE, stderr = FALSE
    )
    status == 0
  }, logical(1))
}
