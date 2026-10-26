#' Install optional tools for `tidygenclust`
#'
#' The core functionality of `tidygenclust` is native: [gt_openadmixture()]
#' and [gt_clumppling()] need neither Python nor conda. A few additional
#' features rely on external tools, which are optional and are only installed
#' on request:
#' - `fastmixture`, a python package used by [gt_fastmixture()];
#' - `ADMIXTURE`, used by [tidypopgen::gt_admixture()].
#'
#' We use `reticulate` to install them in a conda environment,
#' `ctidygenclust`. Additionally, for silicon Macs, `ADMIXTURE` is installed in
#' a separate conda environment `cadmixture86`, as it is only available for OSX
#' as x86 in bioconda.
#'
#' `tgc_tools_install()` creates the conda environment with the chosen
#' `tools`. If the environment already exists, it only adds the tools that are
#' not installed yet (alternatively, use [tgc_fastmixture_install()] and
#' [tgc_admixture_install()] to add a tool to an existing environment).
#' [tgc_tools_version()] reports the tools that are installed.
#'
#' @details For each tool, default to the latest version of these
#'   packages that has been tested to work with `tidygenclust`. It is possible
#'   to provide a more recent github commit for a specific tool, but this might
#'   lead to incompatibilities and errors.
#'
#'   The conda environment is created from a conda yaml file with the
#'   versions of all the packages that have been tested on each operating
#'   system (see `conda_yaml`); `fastmixture` is then installed from its
#'   github repository with `pip`, run in the environment with
#'   `reticulate::conda_run2()`. We have found installation on OSX to be
#'   tricky: if it fails with the tested yaml, you can try `conda_yaml = NULL`
#'   to let conda install the latest available packages.
#'   For OSX, you might also need to install a suitable compiler for openmp
#'   using `brew` in `bash`, setting the correct paths to use it:
#'
#'   `brew install llvm libomp`
#'
#' @param tools a vector with the tools to install, any of "fastmixture" and
#'   "admixture". Defaults to both.
#' @param reset a boolean used to reset the virtual environment. Only set to
#'   TRUE if you have a broken virtual environment that you want to reset.
#' @param fastmixture_hash a string with the commit hash of the `fastmixture`
#'   version to install. Default is the latest tested version.
#' @param conda_yaml A vector of string with the path for the `yaml` files used
#'   for installation (one for `ctidygenclust` and, for osx only, one for
#'   admixture in x86 compatibility mode). The value of each
#'   string can either be "auto" (the default, select the latest available yaml
#'   that has been tested to work on a given operating system), NULL (let conda
#'   attempt to install the latest available packages from the internet), or a
#'   path to a specific conda yaml found in "inst/env_snapshots". For "auto" or
#'   NULL, it is also possible to provide a vector of length one, and the value
#'   will apply to all the conda environments (e.g. `conda_yaml="auto"` is
#'   equivalent to `conda_yaml = c("auto", "auto")` in osx). Tools that are
#'   not in `tools` are left out of the yaml files.
#' @param ci_install a boolean indicating if the installation is being run on
#'   continuous integration (CI) services. Default is FALSE. If TRUE, the
#'   function will look for the conda yaml file in the `inst/env_snapshots`
#'   folder of
#'   the package source directory, rather than in the installed package
#'   directory. This is useful when testing the package on CI services.
#' @returns NULL (invisibly)
#' @export

tgc_tools_install <-
  function(tools = c("fastmixture", "admixture"),
           reset = FALSE,
           fastmixture_hash = "1f03a8980be4872341e5faca9d6f2dfc7ccbfe35",
           conda_yaml = "auto",
           ci_install = FALSE) {
    tgc_check_os()
    tools <- match.arg(tools, several.ok = TRUE)
    is_mac <- tgc_is_mac()
    # the paths of the yaml files (NA to install the latest packages): the
    # first one for ctidygenclust and, on osx, the second for cadmixture86
    conda_yaml <- tgc_resolve_conda_yaml(conda_yaml, ci_install)

    # for osx, check that we have installed the right packages in brew
    if (is_mac && "fastmixture" %in% tools) {
      tgc_check_brew()
    }

    # ctidygenclust holds all the tools, except ADMIXTURE on osx (so it is not
    # needed on osx if only ADMIXTURE is requested)
    envname <- tgc_conda_env()
    needs_env <- "fastmixture" %in% tools || ("admixture" %in% tools && !is_mac)

    # with reset, remove the existing environments first, so that they are
    # created again from scratch below: ctidygenclust only if it is needed
    # (on osx, resetting ADMIXTURE alone must not remove fastmixture), and
    # cadmixture86 only if ADMIXTURE is requested on osx. The tools that were
    # installed in ctidygenclust are installed again, even if not requested,
    # so that resetting an environment does not remove them
    if (reset) {
      if (needs_env && reticulate::condaenv_exists(envname)) {
        in_env <- c("fastmixture", if (!is_mac) "admixture")
        in_env <- in_env[vapply(in_env, tgc_tool_installed, logical(1))]
        also_reset <- setdiff(in_env, tools)
        if (length(also_reset) > 0) {
          message(
            "The conda environment '", envname, "' also contains ",
            paste(also_reset, collapse = " and "), ", which will be ",
            "installed again"
          )
          tools <- c(tools, also_reset)
        }
        reticulate::conda_remove(envname)
      }
      reset_adm86 <- is_mac && "admixture" %in% tools
      if (reset_adm86 && reticulate::condaenv_exists("cadmixture86")) {
        reticulate::conda_remove("cadmixture86")
      }
    }

    # a new environment is created from the yaml without the tools that were
    # not requested; an existing one is kept, and the missing tools are added
    # below
    new_env <- FALSE
    if (needs_env) {
      if (reticulate::condaenv_exists(envname)) {
        message(
          "The conda environment '", envname, "' already exists: only the ",
          "tools that are not installed yet will be added (use ",
          "'reset = TRUE' to recreate it)"
        )
      } else {
        tgc_create_conda_env(envname, tools, conda_yaml[1])
        new_env <- TRUE
      }
    }

    # fastmixture is always (re)installed from github in a new environment,
    # as the conda yaml installs the version from pypi
    if ("fastmixture" %in% tools) {
      if (new_env || !tgc_tool_installed("fastmixture")) {
        tgc_fastmixture_install(fastmixture_hash = fastmixture_hash)
      } else {
        message(
          "fastmixture is already installed; use tgc_fastmixture_install() ",
          "to reinstall it"
        )
      }
    }

    # on linux, ADMIXTURE is usually already installed from the yaml of a new
    # environment; on osx, cadmixture86 is created by tgc_admixture_install()
    if ("admixture" %in% tools) {
      if (!tgc_tool_installed("admixture")) {
        tgc_admixture_install(
          conda_yaml = if (is_mac) conda_yaml[2] else "auto",
          ci_install = ci_install
        )
      } else {
        message("ADMIXTURE is already installed")
      }
    }
    invisible(NULL)
  }

#' Add `fastmixture` to the conda environment of `tidygenclust`
#'
#' `fastmixture` (used by [gt_fastmixture()]) is an optional tool. This
#' function adds it to an existing `ctidygenclust` environment (created with
#' [tgc_tools_install()]), or reinstalls it. See [tgc_tools_install()] for
#' details on installing it on OSX.
#'
#' @inheritParams tgc_tools_install
#' @returns NULL (invisibly)
#' @export

tgc_fastmixture_install <-
  function(fastmixture_hash = "1f03a8980be4872341e5faca9d6f2dfc7ccbfe35") {
    tgc_check_os()
    envname <- tgc_check_conda_env("fastmixture")
    is_mac <- tgc_is_mac()
    if (is_mac) {
      tgc_check_brew()
    }

    # create command line to install fastmixture (from the given commit of
    # its github repository, replacing any version already installed)
    # threadpoolctl is used by gt_fastmixture() to limit the number of
    # threads of fastmixture. --no-deps stops --force-reinstall from also
    # reinstalling the dependencies from pypi (numpy, the only one of
    # fastmixture, is already in the environment, pinned by the conda yaml)
    fast_install_cmd <- paste0(
      "python -m pip install ",
      "--upgrade --force-reinstall --no-deps ",
      "git+https://github.com/Rosemeis/fastmixture.git@",
      fastmixture_hash, " threadpoolctl"
    )

    # on OSX we need to also install a suitable compiler for openmp
    if (is_mac) {
      # install clang and llvm-openmp in the environment, and run pip with
      # the brew llvm first in the PATH and openmp enabled, so that the cython
      # extensions of fastmixture are compiled with openmp
      reticulate::conda_install(
        envname = envname,
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

    # now install fastmixture in the conda environment (conda_run2() writes
    # the lines of cmd_line to a script, run after activating the environment)
    reticulate::conda_run2(
      cmd_line = fast_install_cmd,
      envname = envname
    )
    if (!tgc_tool_installed("fastmixture")) {
      warning(
        "fastmixture has not been successfully installed in the conda ",
        "environment '", envname, "'"
      )
    }
    invisible(NULL)
  }

#' Add `ADMIXTURE` to the conda environment of `tidygenclust`
#'
#' `ADMIXTURE` (used by [tidypopgen::gt_admixture()]) is an optional tool. On
#' Linux, this function adds it to an existing `ctidygenclust` environment
#' (created with [tgc_tools_install()]). On silicon Macs, `ADMIXTURE` is only
#' available as x86 in bioconda, so it is installed in a separate conda
#' environment, `cadmixture86`, which is created by this function.
#'
#' @param conda_yaml only used on OSX: the path of the `yaml` file used to
#'   create the `cadmixture86` environment, "auto" (the default, the latest
#'   tested yaml) or NULL (the latest version of `ADMIXTURE` in bioconda). See
#'   [tgc_tools_install()].
#' @inheritParams tgc_tools_install
#' @returns NULL (invisibly)
#' @export

tgc_admixture_install <- function(conda_yaml = "auto", ci_install = FALSE) {
  tgc_check_os()
  if (tgc_is_mac()) {
    # cadmixture86 is always created from scratch
    if (reticulate::condaenv_exists("cadmixture86")) {
      reticulate::conda_remove("cadmixture86")
    }
    if (is.null(conda_yaml)) {
      conda_yaml <- NA_character_
    }
    # resolve the yaml of cadmixture86, the second environment on osx (the
    # first value, for ctidygenclust, is a placeholder)
    conda_yaml <- tgc_resolve_conda_yaml(c("auto", conda_yaml), ci_install)[2]
    # ADMIXTURE is only available for osx as x86 in bioconda
    # so we have to create a new environment and set it to x86_64
    reticulate::conda_create("cadmixture86",
      channel = c("bioconda", "conda-forge", "defaults")
    )
    # make conda install x86 (osx-64) packages in this environment from now
    # on (they run with Rosetta on silicon Macs)
    reticulate::conda_run2(
      cmd = "conda",
      args = "config --env --set subdir osx-64",
      envname = "cadmixture86"
    )
    if (!is.na(conda_yaml)) {
      # if a conda yaml is provided, use it to install admixture in the new
      # environment (with the conda binary used by reticulate, which may not
      # be on the PATH, e.g. in RStudio)
      system2(
        command = reticulate::conda_binary(),
        args = c(
          "env", "update", "-n", "cadmixture86", "-f", shQuote(conda_yaml)
        )
      )
    } else {
      # install admixture in the new environment
      reticulate::conda_install(
        envname = "cadmixture86",
        packages = c("admixture"),
        channel = c("bioconda")
      )
    }
  } else {
    # on linux we install admixture in ctidygenclust
    envname <- tgc_check_conda_env("admixture")
    reticulate::conda_install(
      envname = envname,
      packages = c("admixture"),
      channel = c("bioconda")
    )
  }
  if (!tgc_tool_installed("admixture")) {
    warning("ADMIXTURE has not been successfully installed")
  }
  invisible(NULL)
}

# the name of the conda environment with the optional tools; it can be changed
# with an option only to test the installation without touching the default
# environment (tidypopgen::gt_admixture() looks for 'ctidygenclust')
tgc_conda_env <- function() {
  getOption("tidygenclust.conda_env", default = "ctidygenclust")
}

# the conda environment of an optional tool: ADMIXTURE has its own
# environment on osx
tgc_tool_env <- function(tool) {
  if (tool == "admixture" && tgc_is_mac()) {
    "cadmixture86"
  } else {
    tgc_conda_env()
  }
}

# check whether an optional tool is installed, by looking for its executable
# in its conda environment (without starting python or conda)
tgc_tool_installed <- function(tool) {
  bin_dir <- tgc_conda_bin(tgc_tool_env(tool))
  !is.null(bin_dir) && file.exists(file.path(bin_dir, tool))
}

# the directory with the executables of a conda environment, or NULL if conda
# or the environment do not exist
tgc_conda_bin <- function(envname) {
  envs <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.null(envs) || !envname %in% envs$name) {
    return(NULL)
  }
  dirname(envs$python[match(envname, envs$name)])
}

# check that the conda environment exists before adding a tool to it, and
# return its name
tgc_check_conda_env <- function(tool) {
  envname <- tgc_conda_env()
  if (!reticulate::condaenv_exists(envname)) {
    stop(
      "The conda environment '", envname, "' does not exist. Create it ",
      "with tgc_tools_install(tools = \"", tool, "\")"
    )
  }
  envname
}

# whether we are on osx (a function, so that tests can mock it)
tgc_is_mac <- function() {
  Sys.info()[["sysname"]] == "Darwin"
}

# give error for windows
tgc_check_os <- function() {
  if (.Platform$OS.type == "windows") {
    stop(
      "tidygenclust does not work on windows; use the Windows subsystem",
      "for Linux (WSL) instead"
    )
  }
}

# on osx, fastmixture needs a compiler for openmp installed with brew
tgc_check_brew <- function() {
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

# find the conda yaml files to use for installation: a vector with the yaml
# for ctidygenclust and, on osx, the one for cadmixture86 (NA if conda should
# install the latest available packages)
tgc_resolve_conda_yaml <- function(conda_yaml, ci_install) {
  is_mac <- tgc_is_mac()
  n_envs <- if (is_mac) 2 else 1
  if (is.null(conda_yaml)) {
    conda_yaml <- NA_character_
  }
  conda_yaml <- as.character(conda_yaml)
  conda_yaml[is.na(conda_yaml) | conda_yaml == "NULL"] <- NA_character_
  # if conda_yaml is a vector of length 1, repeat it to match the number of
  # conda environments
  if (length(conda_yaml) == 1) {
    conda_yaml <- rep(conda_yaml, n_envs)
  }

  # the snapshots are in inst/env_snapshots of the source package, which is
  # in the working directory on CI, or in env_snapshots once installed
  conda_yaml_dir <- ifelse(
    ci_install,
    file.path("inst", "env_snapshots"),
    system.file("env_snapshots/", package = "tidygenclust")
  )
  # replace "auto" with the path to the latest available yaml (the file
  # names end with the date of the snapshot, as YYYYMMDD, so the latest one
  # is the last in alphabetical order)
  patterns <- if (is_mac) {
    c("ctidygenclust_osx_", "cadmixture86_osx_")
  } else {
    "ctidygenclust_linux_"
  }
  for (i in seq_len(n_envs)) {
    if (!is.na(conda_yaml[i]) && conda_yaml[i] == "auto") {
      conda_yaml[i] <- utils::tail(
        list.files(conda_yaml_dir, pattern = patterns[i], full.names = TRUE),
        1
      )
    }
  }

  # check that the files exist (NA values are left as they are)
  for (i in seq_len(n_envs)) {
    if (!is.na(conda_yaml[i]) && !file.exists(conda_yaml[i])) {
      stop(
        paste0(
          "The conda yaml file '",
          conda_yaml[i],
          "' does not exist. Please provide a valid path to a ",
          "conda yaml file."
        )
      )
    }
  }
  conda_yaml[seq_len(n_envs)]
}

# create the conda environment for the optional tools, from a yaml file (with
# the tools that were not requested removed) or from the latest packages
tgc_create_conda_env <- function(envname, tools, conda_yaml) {
  if (!is.na(conda_yaml)) {
    message("Using the conda yaml file ", conda_yaml)
    # use a filtered copy of the yaml, deleted on exit
    conda_yaml <- tgc_filter_conda_yaml(conda_yaml, tools)
    on.exit(unlink(conda_yaml))
  } else {
    conda_yaml <- NULL
  }
  # when a yaml file is given (environment), reticulate creates the
  # environment from it, and 'packages' and 'channel' are ignored
  reticulate::conda_create(
    envname = envname,
    packages = c("python==3.11", "numpy>2.0.0", "cython>3.0.0", "pip"),
    channel = c("bioconda", "conda-forge", "defaults"),
    environment = conda_yaml
  )
}

# a copy of a conda yaml file without the tools that were not requested
tgc_filter_conda_yaml <- function(conda_yaml, tools) {
  lines <- readLines(conda_yaml)
  # the dependency of each tool, either a conda package (e.g.
  # "  - admixture=1.3.0") or a pip one (e.g. "      - fastmixture==1.3.0"),
  # with or without a version
  tool_lines <- c(
    fastmixture = "^\\s*- fastmixture\\s*([=<>]|$)",
    admixture = "^\\s*- admixture\\s*([=<>]|$)"
  )
  for (tool in setdiff(names(tool_lines), tools)) {
    lines <- lines[!grepl(tool_lines[[tool]], lines)]
  }
  # remove the pip section if it is now empty (i.e. it is not followed by
  # more indented lines)
  indent <- function(x) nchar(sub("^(\\s*).*$", "\\1", x))
  pip_line <- grep("^\\s*- pip:\\s*$", lines)
  if (length(pip_line) == 1) {
    next_line <- lines[pip_line + 1]
    if (is.na(next_line) || indent(next_line) <= indent(lines[pip_line])) {
      lines <- lines[-pip_line]
    }
  }
  filtered_yaml <- tempfile(fileext = ".yml")
  writeLines(lines, filtered_yaml)
  filtered_yaml
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
