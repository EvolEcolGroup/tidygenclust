# Tests of the helpers used to install the optional tools (the installation
# itself needs internet access and takes several minutes, so it is not
# tested here)

test_that("conda yaml files are filtered by the requested tools", {
  linux_yaml <- utils::tail(list.files(
    system.file("env_snapshots", package = "tidygenclust"),
    pattern = "ctidygenclust_linux_", full.names = TRUE
  ), 1)
  yaml_lines <- readLines(linux_yaml)
  has_line <- function(lines, pattern) any(grepl(pattern, lines))
  # all the tools: the file is unchanged
  filtered <- tgc_filter_conda_yaml(linux_yaml, c("fastmixture", "admixture"))
  on.exit(unlink(filtered), add = TRUE)
  expect_identical(readLines(filtered), yaml_lines)
  # only fastmixture
  filtered <- tgc_filter_conda_yaml(linux_yaml, "fastmixture")
  on.exit(unlink(filtered), add = TRUE)
  filtered_lines <- readLines(filtered)
  expect_false(has_line(filtered_lines, "- admixture="))
  expect_true(has_line(filtered_lines, "- fastmixture=="))
  expect_true(has_line(filtered_lines, "- pip:"))
  expect_length(filtered_lines, length(yaml_lines) - 1)
  # only admixture: the pip section, which only had fastmixture, is removed
  filtered <- tgc_filter_conda_yaml(linux_yaml, "admixture")
  on.exit(unlink(filtered), add = TRUE)
  filtered_lines <- readLines(filtered)
  expect_true(has_line(filtered_lines, "- admixture="))
  expect_false(has_line(filtered_lines, "fastmixture"))
  expect_false(has_line(filtered_lines, "- pip:"))
  expect_true(has_line(filtered_lines, "- python="))
  expect_length(filtered_lines, length(yaml_lines) - 2)
  # a pip section with other packages is kept
  yaml_file <- withr::local_tempfile(fileext = ".yml")
  writeLines(c(
    "name: test", "dependencies:", "  - pip", "  - pip:",
    "      - fastmixture==1.3.0", "      - other==1.0", "prefix: /x"
  ), yaml_file)
  filtered <- tgc_filter_conda_yaml(yaml_file, "admixture")
  on.exit(unlink(filtered), add = TRUE)
  expect_identical(readLines(filtered), c(
    "name: test", "dependencies:", "  - pip", "  - pip:",
    "      - other==1.0", "prefix: /x"
  ))
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
# environments and the installed tools (as "tool@environment").
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
    tgc_check_brew = function() NULL,
    tgc_conda_env = function() "ctidygenclust",
    tgc_tool_installed = function(tool) {
      paste0(tool, "@", tool_env(tool)) %in% state$installed
    },
    tgc_create_conda_env = function(envname, tools, conda_yaml) {
      record(paste("create", envname, paste(sort(tools), collapse = "+")))
      state$envs <- c(state$envs, envname)
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
  expect_message(tgc_tools_install(), "already exists")
  expect_identical(state$log, "install admixture")
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
