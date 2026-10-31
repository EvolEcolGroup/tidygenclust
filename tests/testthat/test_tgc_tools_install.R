# Tests of the helpers used to install the optional tools (the installation
# itself needs internet access and takes several minutes, so it is not
# tested here)

test_that("conda yaml files are filtered by the requested tools", {
  # a snapshot in the format of those in inst/env_snapshots, written here so
  # that the test does not depend on the content of the latest snapshot
  yaml_lines <- c(
    "name: ctidygenclust",
    "channels:",
    "  - conda-forge",
    "  - bioconda",
    "dependencies:",
    "  - admixture=1.3.0",
    "  - numpy=2.3.3",
    "  - pip=26.2.1",
    "  - python=3.11.0",
    "  - pip:",
    "      - fastmixture==1.3.0",
    "prefix: /usr/share/miniconda/envs/ctidygenclust"
  )
  linux_yaml <- withr::local_tempfile(fileext = ".yml")
  writeLines(yaml_lines, linux_yaml)
  filter_lines <- function(yaml, tools) {
    filtered <- tgc_filter_conda_yaml(yaml, tools)
    on.exit(unlink(filtered))
    readLines(filtered)
  }
  # all the tools: the file is unchanged
  expect_identical(
    filter_lines(linux_yaml, c("fastmixture", "admixture")), yaml_lines
  )
  # only fastmixture
  expect_identical(filter_lines(linux_yaml, "fastmixture"), yaml_lines[-6])
  # only admixture: the pip section, which only had fastmixture, is removed
  expect_identical(filter_lines(linux_yaml, "admixture"), yaml_lines[-10:-11])
  # a pip section with other packages is kept
  yaml_file <- withr::local_tempfile(fileext = ".yml")
  writeLines(c(
    "name: test", "dependencies:", "  - pip", "  - pip:",
    "      - fastmixture==1.3.0", "      - other==1.0", "prefix: /x"
  ), yaml_file)
  expect_identical(filter_lines(yaml_file, "admixture"), c(
    "name: test", "dependencies:", "  - pip", "  - pip:",
    "      - other==1.0", "prefix: /x"
  ))
})

test_that("the tools are found in the latest snapshots", {
  # the latest snapshot of ctidygenclust for each operating system: the
  # filter must find the lines of the tools in its format (the snapshots of
  # macOS do not have admixture, which has its own environment)
  snapshots <- list.files(
    system.file("env_snapshots", package = "tidygenclust"),
    pattern = "^ctidygenclust_", full.names = TRUE
  )
  os <- sub("^ctidygenclust_(.*)_\\d+\\.yml$", "\\1", basename(snapshots))
  for (yaml in tapply(snapshots, os, function(x) utils::tail(sort(x), 1))) {
    lines <- readLines(yaml)
    expect_true(any(grepl("^\\s*- fastmixture==", lines)), info = yaml)
    filtered <- tgc_filter_conda_yaml(yaml, character(0))
    on.exit(unlink(filtered), add = TRUE)
    filtered_lines <- readLines(filtered)
    expect_false(any(grepl("fastmixture|- admixture", filtered_lines)),
      info = yaml
    )
    # nothing else is removed, apart from an empty pip section
    removed <- setdiff(lines, filtered_lines)
    expect_true(all(grepl("fastmixture|- admixture|- pip:", removed)),
      info = yaml
    )
  }
})

test_that("conda yaml files are found", {
  n_envs <- if (Sys.info()[["sysname"]] == "Darwin") 2 else 1
  yaml_files <- tgc_resolve_conda_yaml("auto", ci_install = FALSE)
  expect_length(yaml_files, n_envs)
  expect_true(all(file.exists(yaml_files)))
  expect_match(yaml_files[1], "ctidygenclust_")
  expect_identical(
    tgc_resolve_conda_yaml(NULL, ci_install = FALSE),
    rep(NA_character_, n_envs)
  )
  expect_identical(
    tgc_resolve_conda_yaml(yaml_files, ci_install = FALSE),
    yaml_files
  )
  expect_error(
    tgc_resolve_conda_yaml("not_a_file.yml", ci_install = FALSE),
    "does not exist"
  )
})

test_that("missing tools are reported", {
  # the installers stop on windows, before looking for the environment (as
  # tested below), so these tests are of the other systems
  testthat::local_mocked_bindings(tgc_is_windows = function() FALSE)
  withr::local_options(tidygenclust.conda_env = "tidygenclust_no_such_env")
  expect_false(tgc_tool_installed("fastmixture"))
  if (Sys.info()[["sysname"]] != "Darwin") {
    expect_false(tgc_tool_installed("admixture"))
  }
  expect_identical(tgc_tools_version()$fastmixture, NA)
  expect_error(
    gt_fastmixture("plink_prefix", k = 2),
    "tgc_tools_install\\(tools = \"fastmixture\"\\)"
  )
  expect_error(
    tgc_fastmixture_install(),
    "'tidygenclust_no_such_env' does not exist"
  )
  if (Sys.info()[["sysname"]] != "Darwin") {
    expect_error(
      tgc_admixture_install(),
      "'tidygenclust_no_such_env' does not exist"
    )
  }
  expect_error(tgc_tools_install(tools = "clumppling"), "should be one of")
})

# Replace the functions that query and change conda environments with mocks
# that record what the installer does. `state` holds the existing
# environments, the installed tools (as "tool@environment") and the yaml
# files used to create environments.
local_mock_conda <- function(envs, installed, mac = FALSE,
                             env = parent.frame()) {
  state <- new.env()
  state$envs <- envs
  state$installed <- installed
  state$log <- character()
  record <- function(x) state$log <- c(state$log, x)
  tool_env <- function(tool) {
    if (tool == "admixture" && mac) "cadmixture86" else "ctidygenclust"
  }
  install <- function(tool) {
    record(paste("install", tool))
    state$installed <- c(state$installed, paste0(tool, "@", tool_env(tool)))
  }
  testthat::local_mocked_bindings(
    condaenv_exists = function(envname = NULL, ...) envname %in% state$envs,
    conda_remove = function(envname, ...) {
      record(paste("remove", envname))
      state$envs <- setdiff(state$envs, envname)
      state$installed <- state$installed[
        !endsWith(state$installed, paste0("@", envname))
      ]
    },
    .package = "reticulate", .env = env
  )
  testthat::local_mocked_bindings(
    tgc_is_mac = function() mac,
    tgc_is_windows = function() FALSE,
    tgc_check_brew = function() NULL,
    tgc_conda_env = function() "ctidygenclust",
    tgc_tool_installed = function(tool) {
      paste0(tool, "@", tool_env(tool)) %in% state$installed
    },
    tgc_create_conda_env = function(envname, tools, conda_yaml) {
      record(paste("create", envname, paste(sort(tools), collapse = "+")))
      state$envs <- c(state$envs, envname)
      state$yaml <- c(state$yaml, conda_yaml)
    },
    tgc_fastmixture_install = function(...) install("fastmixture"),
    tgc_admixture_install = function(...) install("admixture"),
    .env = env
  )
  state
}

test_that("resetting the environment keeps the other tools", {
  # on linux, both tools are in ctidygenclust: resetting it to reinstall
  # ADMIXTURE also installs fastmixture again
  state <- local_mock_conda(
    envs = "ctidygenclust",
    installed = c("fastmixture@ctidygenclust", "admixture@ctidygenclust")
  )
  expect_message(
    tgc_tools_install(tools = "admixture", reset = TRUE),
    "also contains fastmixture"
  )
  expect_identical(state$log, c(
    "remove ctidygenclust", "create ctidygenclust admixture+fastmixture",
    "install fastmixture", "install admixture"
  ))
  # on osx, ADMIXTURE is in its own environment: ctidygenclust (with
  # fastmixture) must be left alone
  state <- local_mock_conda(
    envs = c("ctidygenclust", "cadmixture86"),
    installed = c("fastmixture@ctidygenclust", "admixture@cadmixture86"),
    mac = TRUE
  )
  tgc_tools_install(tools = "admixture", reset = TRUE)
  expect_identical(state$log, c("remove cadmixture86", "install admixture"))
  # without reset, only the missing tools are added
  state <- local_mock_conda(
    envs = "ctidygenclust",
    installed = "fastmixture@ctidygenclust"
  )
  msgs <- capture_messages(tgc_tools_install())
  expect_match(msgs, "already exists", all = FALSE)
  expect_identical(state$log, "install admixture")
})

test_that("new environments are created for the requested tools only", {
  # linux, fastmixture only: ctidygenclust is created from the latest
  # snapshot (filtered by tgc_create_conda_env()), and fastmixture is
  # installed from github
  state <- local_mock_conda(envs = character(), installed = character())
  tgc_tools_install(tools = "fastmixture")
  expect_identical(state$log, c(
    "create ctidygenclust fastmixture", "install fastmixture"
  ))
  expect_match(basename(state$yaml), "^ctidygenclust_linux_\\d+\\.yml$")
  # without a yaml file, the latest packages are installed
  state <- local_mock_conda(envs = character(), installed = character())
  tgc_tools_install(tools = "fastmixture", conda_yaml = NULL)
  expect_identical(state$yaml, NA_character_)
  # osx, ADMIXTURE only: ctidygenclust is not needed
  state <- local_mock_conda(
    envs = character(), installed = character(), mac = TRUE
  )
  tgc_tools_install(tools = "admixture")
  expect_identical(state$log, "install admixture")
  # osx, both tools (the yaml of ctidygenclust for osx has no ADMIXTURE, so
  # filtering it for ADMIXTURE keeps nothing more)
  state <- local_mock_conda(
    envs = character(), installed = character(), mac = TRUE
  )
  tgc_tools_install()
  expect_identical(state$log, c(
    "create ctidygenclust admixture+fastmixture", "install fastmixture",
    "install admixture"
  ))
  expect_match(basename(state$yaml), "^ctidygenclust_osx_\\d+\\.yml$")
  # everything installed already: nothing is done
  state <- local_mock_conda(
    envs = "ctidygenclust",
    installed = c("fastmixture@ctidygenclust", "admixture@ctidygenclust")
  )
  msgs <- capture_messages(tgc_tools_install())
  expect_match(msgs, "'ctidygenclust' already exists", all = FALSE)
  expect_match(msgs, "fastmixture is already installed", all = FALSE)
  expect_match(msgs, "ADMIXTURE is already installed", all = FALSE)
  expect_identical(state$log, character())
})

test_that("fastmixture is installed without reinstalling its dependencies", {
  # --force-reinstall alone also reinstalled numpy from pypi, replacing the
  # version pinned by the conda yaml
  run <- new.env()
  testthat::local_mocked_bindings(
    conda_run2 = function(cmd_line, envname, ...) run$cmd_line <- cmd_line,
    .package = "reticulate"
  )
  testthat::local_mocked_bindings(
    tgc_is_mac = function() FALSE,
    tgc_is_windows = function() FALSE,
    tgc_check_conda_env = function(tool) "ctidygenclust",
    tgc_tool_installed = function(tool) TRUE
  )
  tgc_fastmixture_install(fastmixture_hash = "abc123")
  expect_length(run$cmd_line, 1)
  expect_match(run$cmd_line, "pip install .*--force-reinstall --no-deps")
  expect_match(run$cmd_line, "fastmixture.git@abc123 threadpoolctl$")
})

test_that("on osx, ADMIXTURE is installed with reticulate's conda", {
  # fake conda executables, which write their arguments to a file: one is
  # the binary used by reticulate, the other is first on the PATH (the old
  # code called the conda on the PATH, which may not exist, e.g. in RStudio)
  fake_dir <- withr::local_tempdir()
  fake_conda <- function(name, out) {
    path <- file.path(fake_dir, name)
    writeLines(c("#!/bin/sh", paste0("printf '%s\\n' \"$@\" > '", out, "'")),
      path
    )
    Sys.chmod(path, "755")
    path
  }
  out_reticulate <- file.path(fake_dir, "args_reticulate.txt")
  out_path <- file.path(fake_dir, "args_path.txt")
  reticulate_conda <- fake_conda("reticulate_conda", out_reticulate)
  dir.create(file.path(fake_dir, "bin"))
  fake_conda(file.path("bin", "conda"), out_path)
  withr::local_path(file.path(fake_dir, "bin"), action = "prefix")
  # a yaml file whose path contains a space
  yaml_file <- file.path(fake_dir, "my env.yml")
  writeLines("name: cadmixture86", yaml_file)
  run <- new.env()
  testthat::local_mocked_bindings(
    condaenv_exists = function(envname = NULL, ...) FALSE,
    conda_create = function(envname, ...) NULL,
    conda_run2 = function(cmd, args, envname, ...) run$args <- args,
    conda_binary = function(...) reticulate_conda,
    .package = "reticulate"
  )
  testthat::local_mocked_bindings(
    tgc_is_mac = function() TRUE,
    tgc_is_windows = function() FALSE,
    tgc_tool_installed = function(tool) TRUE
  )
  tgc_admixture_install(conda_yaml = yaml_file)
  expect_identical(run$args, "config --env --set subdir osx-64")
  expect_false(file.exists(out_path))
  expect_identical(
    readLines(out_reticulate),
    c("env", "update", "-n", "cadmixture86", "-f", yaml_file)
  )
})

test_that("on osx, a single yaml file is only used for ctidygenclust", {
  testthat::local_mocked_bindings(tgc_is_mac = function() TRUE)
  my_yaml <- withr::local_tempfile(fileext = ".yml")
  writeLines("name: ctidygenclust", my_yaml)
  # the yaml of ctidygenclust must not be used for cadmixture86
  yaml_files <- tgc_resolve_conda_yaml(my_yaml, ci_install = FALSE)
  expect_identical(yaml_files[1], my_yaml)
  expect_match(yaml_files[2], "cadmixture86_osx_")
  # "auto" and NULL apply to both environments
  yaml_files <- tgc_resolve_conda_yaml("auto", ci_install = FALSE)
  expect_match(yaml_files[1], "ctidygenclust_osx_")
  expect_match(yaml_files[2], "cadmixture86_osx_")
  expect_identical(
    tgc_resolve_conda_yaml(NULL, ci_install = FALSE),
    c(NA_character_, NA_character_)
  )
})

test_that("windows is rejected with a readable message", {
  testthat::local_mocked_bindings(tgc_is_windows = function() TRUE)
  expect_error(
    tgc_tools_install(),
    "use the Windows subsystem for Linux \\(WSL\\) instead"
  )
  # the other installers are rejected too, saying that only the optional
  # tools are not available
  expect_error(tgc_fastmixture_install(), "cannot be installed on windows")
  expect_error(tgc_admixture_install(), "cannot be installed on windows")
})

test_that("a missing snapshot directory gives a clear error", {
  # with ci_install, the snapshots are looked for in inst/env_snapshots of
  # the working directory
  withr::local_dir(withr::local_tempdir())
  expect_error(
    tgc_resolve_conda_yaml("auto", ci_install = TRUE),
    "No conda yaml snapshot found .*run from the root of the package source"
  )
})

test_that("on osx, fastmixture is compiled with the llvm of homebrew", {
  run <- new.env()
  testthat::local_mocked_bindings(
    conda_install = function(...) NULL,
    conda_run2 = function(cmd_line, envname, ...) run$cmd_line <- cmd_line,
    .package = "reticulate"
  )
  testthat::local_mocked_bindings(
    tgc_is_mac = function() TRUE,
    tgc_is_windows = function() FALSE,
    tgc_check_brew = function() NULL,
    tgc_check_conda_env = function(tool) "ctidygenclust",
    tgc_tool_installed = function(tool) TRUE
  )
  tgc_fastmixture_install()
  # the location of llvm is asked to homebrew (it was hard-coded for silicon
  # Macs, /opt/homebrew)
  expect_identical(
    run$cmd_line[1],
    "export PATH=\"$(brew --prefix llvm)/bin:$PATH\""
  )
  expect_match(run$cmd_line[length(run$cmd_line)], "pip install")
})
