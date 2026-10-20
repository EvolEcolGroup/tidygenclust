#' Install tools for `tidygenclust`
#'
#' `tidygenclust` relies on `ADMIXTURE` and on python packages `fastmixture`
#' and `clumppling` for a number of
#' functionalities. On Linux and macOS we use `reticulate` to install them in
#' conda environments; on Windows we use micromamba inside WSL2 for the Python
#' tools only.
#' As their dependencies are incompatible, we use two separate
#' conda environments, `ctidygenclust` (for `fastmixture` and `admixture`) and
#' `cclumppling` (for `clumppling`). Additionally, for silicon Macs, `ADMIXTURE`
#' is installed in a separate conda environment `cadmixture86`, as it is only
#' available for OSX as x86 in bioconda.
#' @details
#' For each tool, default to the latest tested version of
#' these packages that have been tested to work with `tidegenclust`. It is
#' possible to provide a more recent github commit for a specific tool, but
#' this might lead to incompatibilities and errors.
#'
#' We have found installation on OSX to be tricky, so we provide two methods
#' for installing `fastmixture` on OSX: `reticulate` and `conda_yaml`. The
#' `reticulate` method uses the `reticulate::conda_run2()` function to run
#' installation commands, while the `conda_yaml` method creates a conda
#' environment directly with conda. If the `reticulate` method fails, you can
#' use the `conda_yaml` method to create the environment directly with conda.
#' For OSX, you might also need to install a suitable compiler for openmp
#' using `brew` in `bash`, setting the correct paths to use it:
#'
#' `brew install llvm libomp`
#'
#' `ADAMIXTURE` (used by [gt_adamixture()]) is optional. Set
#' `adamixture = TRUE` to install it
#' in `ctidygenclust` together with the other tools, or use
#' [tgc_adamixture_install()] to add it to an existing `ctidygenclust`
#' environment.
#'
#' Windows support is experimental and requires an x86-64 computer with WSL2
#' and Ubuntu. The tools are installed in user-level environments inside WSL2
#' and invoked automatically from the Windows R session. The default WSL
#' distribution is `Ubuntu`; set
#' `options(tidygenclust.wsl_distro = "DistributionName")` before installation
#' to use another distribution. Alternatively, set `TIDYGENCLUST_WSL_DISTRO`
#' in your Windows `.Renviron` file; this also applies to new R processes.
#' The R option takes precedence over the environment variable.
#' Re-running the Windows installer reuses the environments and updates tools
#' whose requested source revisions differ from those installed.
#' Fresh installations and explicit rebuilds use bundled linux-64 dependency
#' locks: exact Conda builds and SHA-256-checked Python dependencies, including
#' source-build tools. Existing environments are not silently rebuilt to match
#' the locks. Custom source revisions must work with these dependencies.
#' This Windows backend supports Fastmixture, ADAMIXTURE and Clumppling only;
#' it does not extend [tidypopgen::gt_admixture()] to Windows.
#' See `vignette("windows", package = "tidygenclust")` for setup instructions.
#'
#' @param reset a boolean used to reset the virtual environment. Only set
#' to TRUE if you have a broken virtual environment that you want to reset.
#' @param fastmixture_hash a string with the commit hash of the `fastmixture`
#' version to install. Default is the latest tested version.
#' @param clumppling_hash a string with the commit hash of the `clumppling`
#' version to install. Default is the latest tested version.
#' @param conda_method a string indicating the method to create the environment
#' used for `fastmixture`. Default is `reticulate`, which uses the
#' `reticulate::conda_run2()` function to run the installation commands (this is
#' the default, and the only method for Linux.
#' Alternatively, for OSX, you can use `conda_yaml`, which will create a conda
#' environment directly with conda. Use this second method if "reticulate" fails
#' whilst trying to install on OSX.
#' @param ci_install a boolean indicating if the installation is being run on
#' continuous integration (CI) services. Default is FALSE. If TRUE, the
#' function will look for the conda yaml file in the `inst/python` folder
#' of the package source directory, rather than in the installed package
#' directory. This is useful when testing the package on CI services.
#' @param adamixture a boolean indicating whether to install `ADAMIXTURE` in
#' `ctidygenclust`. Default is FALSE.
#' @param adamixture_hash a string with the commit hash of the `ADAMIXTURE`
#' version to install. Default is the latest tested version.
#' @returns NULL
#' @export

tgc_tools_install <-
  function(reset = FALSE,
           fastmixture_hash = "29e04339ce6ddf750ee4e06f8aabe40335e0d0ee",
           clumppling_hash = "2d24e0b2f6ddfcb51a436df96a06d5f57d18d20a",
           conda_method = c("reticulate", "conda_yaml"),
           ci_install = FALSE,
           adamixture = FALSE,
           adamixture_hash = "43dca994a6091ae7034fbec900ce7bb9a5c6d8c1") {
    # Windows installs inside WSL; leave the native Conda installer unchanged.
    if (.Platform$OS.type == "windows") {
      return(.tgc_wsl_install_tools(
        reset = reset,
        fastmixture_hash = fastmixture_hash,
        clumppling_hash = clumppling_hash,
        adamixture = adamixture,
        adamixture_hash = adamixture_hash
      ))
    }

    if (.Platform$OS.type != "unix") {
      stop(
        "tidygenclust supports Linux, macOS, and Windows through WSL2.",
        call. = FALSE
      )
    }

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
    # check that the osx_method is valid
    conda_method <- match.arg(conda_method)
    # if linux, force method to be reticulate
    if (Sys.info()["sysname"] == "Linux") {
      conda_method <- "reticulate"
    }

    # install fastmixture

    # if method is reticulate
    if (conda_method == "reticulate") {
      # create a conda environment with the necessary packages
      reticulate::conda_create(
        envname = "ctidygenclust",
        packages = c("python==3.11", "numpy>2.0.0", "cython>3.0.0", "pip"),
        channel = c("bioconda", "conda-forge", "defaults")
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

      ## https://github.com/rstudio/reticulate/issues/905
      reticulate::conda_run2(
        cmd_line = fast_install_cmd,
        envname = "ctidygenclust"
      )
    } else if (conda_method == "conda_yaml") {
      # create a conda environment with the necessary packages
      # using a conda yml file
      if (!ci_install) {
        yml_path <- system.file("python/env_osx.yml", package = "tidygenclust")
      } else {
        yml_path <- "inst/python/env_osx.yml"
      }

      cat("yaml path is :\n")
      cat(yml_path)
      cat("\n")

      system(
        command = paste("export PATH=\"/opt/homebrew/opt/llvm/bin:$PATH\"; conda env create -f '", # nolint
          yml_path, "'",
          collapse = "", sep = ""
        )
      )

      # system2(command = "conda",
      #   args = paste("env create -f '",
      #                yml_path,
      #                "'", collapse = "", sep = ""
      #   )
      # )
    }

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


    # optionally install adamixture in ctidygenclust
    if (adamixture) {
      tgc_adamixture_install(adamixture_hash = adamixture_hash)
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
      # Install clumppling
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
        channel = c("bioconda", "conda-forge", "defaults")
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

#' Install `ADAMIXTURE` in the `ctidygenclust` conda environment
#'
#' `ADAMIXTURE` (used by [gt_adamixture()]) is an optional tool.
#' This function adds it to an existing `ctidygenclust` environment
#' (created with [tgc_tools_install()]). It installs a lightweight
#' CPU version of `ADAMIXTURE` without requiring PyTorch.
#' On Windows, dependencies and source-build tools use the same bundled locks
#' as [tgc_tools_install()].
#'
#' @param adamixture_hash a string with the commit hash of the `ADAMIXTURE`
#' version to install. Default is the latest tested version.
#' @param adamixture_repo a string with the github repository (`owner/repo`)
#' from which to install `ADAMIXTURE`.
#' @returns NULL
#' @export

tgc_adamixture_install <-
  function(adamixture_hash = "43dca994a6091ae7034fbec900ce7bb9a5c6d8c1",
           adamixture_repo = "AI-sandbox/ADAMIXTURE") {
    if (.Platform$OS.type == "windows") {
      # Add ADAMIXTURE to the existing core environment, not a new environment.
      if (!.tgc_wsl_env_exists("ctidygenclust")) {
        stop(
          "The WSL environment 'ctidygenclust' does not exist. Create it ",
          "first with 'tgc_tools_install()'.",
          call. = FALSE
        )
      }
      package_spec <- paste0(
        "git+https://github.com/", adamixture_repo, ".git@", adamixture_hash
      )
      # Use the same locked runtime/build dependencies as the main installer.
      # Disable dependency resolution and build isolation to preserve those 
      # pins.
      dependencies <- system.file(
        "wsl", "ctidygenclust-pip.lock",
        package = "tidygenclust", mustWork = TRUE
      )
      .tgc_wsl_run_env(
        "ctidygenclust", "python",
        c(
          "-m", "pip", "install", "--no-cache-dir", "--require-hashes",
          "--no-deps",
          "--no-build-isolation", "--requirement", .tgc_wsl_path(dependencies)
        )
      )
      .tgc_wsl_run_env(
        "ctidygenclust",
        "python",
        c(
          "-m", "pip", "install", "--no-cache-dir", "--force-reinstall",
          "--no-deps",
          "--no-build-isolation", package_spec
        )
      )
      .tgc_wsl_run_env("ctidygenclust", "python", c("-m", "pip", "check"))
      if (!.tgc_tool_available("adamixture")) {
        warning(
          "ADAMIXTURE has not been successfully installed in WSL.",
          call. = FALSE
        )
      }
      return(invisible(NULL))
    }

    if (!reticulate::condaenv_exists("ctidygenclust")) {
      stop(
        "The conda environment 'ctidygenclust' does not exist. Create it ",
        "first with 'tgc_tools_install()'"
      )
    }
    ada_install_cmd <- paste0(
      "python -m pip install ",
      "git+https://github.com/", adamixture_repo, ".git@",
      adamixture_hash
    )
    if (Sys.info()["sysname"] == "Darwin") {
      ada_install_cmd <- c(
        "export PATH=\"/opt/homebrew/opt/llvm/bin:$PATH\"",
        "export CC=clang",
        "export CXX=clang++",
        ada_install_cmd
      )
    }
    out <- reticulate::conda_run2(
      cmd_line = ada_install_cmd,
      envname = "ctidygenclust"
    )
    if (out != 0) {
      warning(
        "ADAMIXTURE has not been successfully installed ",
        "in your conda environment"
      )
    }
    invisible(NULL)
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
