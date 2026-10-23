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
