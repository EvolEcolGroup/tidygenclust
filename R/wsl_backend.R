# Internal helpers for running the analysis tools in WSL from Windows R.
# R keeps the analysis workflow; separate Linux processes run the Python tools.

# An R option overrides the persistent .Renviron setting; Ubuntu is the default.
.tgc_wsl_distro <- function() {
  getOption(
    "tidygenclust.wsl_distro",
    Sys.getenv("TIDYGENCLUST_WSL_DISTRO", "Ubuntu")
  )
}

# Use one entry point for WSL commands, with the usual system2 output controls.
.tgc_wsl_call <- function(args, stdout = "", stderr = "") {
  if (!nzchar(Sys.which("wsl.exe"))) {
    stop(
      "WSL is not installed or 'wsl.exe' is not available on PATH.",
      call. = FALSE
    )
  }

  # system2 constructs a Windows command line: protect spaces and
  # metacharacters.
  # --exec passes the resulting arguments to Linux without a Linux shell.
  args <- as.character(c(
    "--distribution", .tgc_wsl_distro(), "--exec", args
  ))
  needs_quotes <- !nzchar(args) | grepl("[[:space:]&|<>^()\"]", args)
  args[needs_quotes] <- vapply(
    args[needs_quotes], shQuote, character(1),
    type = "cmd"
  )
  # Callers inspect the exit status and supply a more useful error message.
  suppressWarnings(system2(
    command = "wsl.exe",
    args = args,
    stdout = stdout,
    stderr = stderr
  ))
}

# Capture diagnostic output and normalise system2's success status to zero.
.tgc_wsl_capture <- function(args) {
  out <- .tgc_wsl_call(args, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = as.integer(status), output = out)
}

# Check that the selected distribution starts, not just that wsl.exe exists.
# The installer separately checks WSL2 and the supported processor architecture.
.tgc_wsl_available <- function() {
  if (.Platform$OS.type != "windows" || !nzchar(Sys.which("wsl.exe"))) {
    return(FALSE)
  }
  check <- tryCatch(
    .tgc_wsl_capture(c("/usr/bin/test", "-d", "/proc")),
    error = function(e) list(status = 1L)
  )
  identical(check$status, 0L)
}

# Keep the original reticulate route on Linux/macOS; never use Windows Python.
.tgc_backend <- function() {
  if (.Platform$OS.type == "windows") {
    if (.tgc_wsl_available()) {
      return("wsl")
    }
    return("unavailable")
  }
  "reticulate"
}

.tgc_wsl_root <- function() {
  # Tests can use a disposable root without changing the user's environments.
  root <- getOption("tidygenclust.wsl_root")
  if (!is.null(root)) {
    return(sub("/+$", "", root))
  }
  # Read HOME from the selected distribution, not the Windows R session.
  result <- .tgc_wsl_capture(c("/usr/bin/printenv", "HOME"))
  if (result$status != 0L || length(result$output) == 0L) {
    stop(
      "Could not determine the home directory in WSL distribution '",
      .tgc_wsl_distro(), "'.",
      call. = FALSE
    )
  }
  paste0(trimws(result$output[[1]]), "/.local/share/tidygenclust")
}

# Look only in the package's own prefix, not unrelated user Conda installations.
.tgc_wsl_env_exists <- function(envname, root = .tgc_wsl_root()) {
  check <- tryCatch(
    .tgc_wsl_capture(c(
      "/usr/bin/test", "-x", paste0(root, "/envs/", envname, "/bin/python")
    )),
    error = function(e) list(status = 1L)
  )
  identical(check$status, 0L)
}

# Run in an explicit Linux environment without activating it in a shell.
.tgc_wsl_run_env <- function(
  envname,
  command,
  args = character(),
  stdout = "",
  stderr = ""
) {
  # Resolve the paths once: each lookup of the Linux HOME starts wsl.exe.
  root <- .tgc_wsl_root()
  if (!.tgc_wsl_env_exists(envname, root)) {
    stop(
      "The WSL environment '", envname,
      "' does not exist in '", .tgc_wsl_distro(), "'. ",
      "Install it with 'tgc_tools_install()'.",
      call. = FALSE
    )
  }

  status <- .tgc_wsl_call(
    c(
      paste0(root, "/bin/micromamba"),
      "--root-prefix", paste0(root, "/micromamba"),
      "run", "--prefix", paste0(root, "/envs/", envname),
      command, args
    ),
    stdout = stdout,
    stderr = stderr
  )

  if (identical(stdout, TRUE)) {
    # Captured probes inspect the status attribute themselves (e.g. versions).
    return(status)
  }

  if (!identical(as.integer(status), 0L)) {
    stop(
      "The WSL command failed with exit status ", status, ".",
      call. = FALSE
    )
  }
  invisible(status)
}

# Convert Windows paths for Linux, including output paths not yet created.
.tgc_wsl_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Expected one non-empty Windows path.", call. = FALSE)
  }

  windows_path <- normalizePath(
    path.expand(path),
    winslash = "/",
    mustWork = FALSE
  )
  result <- .tgc_wsl_capture(c(
    "/usr/bin/wslpath", "-a", "-u", windows_path
  ))
  if (result$status != 0L || length(result$output) == 0L) {
    stop(
      "Could not translate the Windows path for WSL: ", windows_path,
      if (length(result$output)) {
        paste0("\n", paste(result$output, collapse = "\n"))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  trimws(result$output[[1]])
}

# Test imports in the matching environment without initialising reticulate.
# Preserve the existing availability checks on Linux and macOS.
.tgc_tool_available <- function(tool) {
  if (.Platform$OS.type == "windows") {
    if (!tool %in% c("fastmixture", "adamixture", "clumppling")) {
      return(FALSE)
    }
    envname <- if (tool == "clumppling") "cclumppling" else "ctidygenclust"
    status <- tryCatch(
      .tgc_wsl_run_env(
        envname, "python", c("-c", paste0("import ", tool)),
        stdout = FALSE, stderr = FALSE
      ),
      error = function(e) 1L
    )
    identical(as.integer(status), 0L)
  } else {
    switch(tool,
      fastmixture = reticulate::condaenv_exists("ctidygenclust") &&
        reticulate::py_module_available("fastmixture"),
      adamixture = reticulate::condaenv_exists("ctidygenclust") &&
        reticulate::py_module_available("adamixture"),
      clumppling = reticulate::condaenv_exists("cclumppling"),
      FALSE
    )
  }
}

# Match the Python runner's float64 format. byrow restores NumPy's row order,
# since R would otherwise fill the matrix column by column.
.tgc_read_binary_matrix <- function(path, shape) {
  shape <- as.integer(shape)
  if (length(shape) != 2L || anyNA(shape) || any(shape < 0L)) {
    stop("The WSL runner returned an invalid matrix shape.", call. = FALSE)
  }

  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  values <- readBin(
    connection,
    what = double(),
    n = prod(shape),
    size = 8L,
    endian = "little"
  )
  if (length(values) != prod(shape)) {
    stop("The WSL runner returned a truncated matrix file.", call. = FALSE)
  }
  matrix(values, nrow = shape[[1]], ncol = shape[[2]], byrow = TRUE)
}

# Exchange files with the runner, then return the same objects as reticulate.
.tgc_wsl_model_run <- function(algorithm, args) {
  # Each run has its own directory, removed on success or error.
  job_dir <- tempfile(pattern = "tidygenclust-wsl-")
  if (!dir.create(job_dir, recursive = TRUE)) {
    stop("Could not create a temporary WSL job directory.", call. = FALSE)
  }
  on.exit(unlink(job_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Only file arguments need translation; keep analysis options unchanged.
  args$bfile <- .tgc_wsl_path(args$bfile)
  if (!is.null(args$supervised)) {
    args$supervised <- .tgc_wsl_path(args$supervised)
  }

  config_path <- file.path(job_dir, "job.json")
  result_path <- file.path(job_dir, "result.json")
  # Exchange small arguments as JSON and large matrices as binary files.
  jsonlite::write_json(
    list(
      algorithm = algorithm,
      output_dir = .tgc_wsl_path(job_dir),
      args = args
    ),
    path = config_path,
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  )

  # The runner reuses the package's existing Python analysis wrappers.
  runner <- system.file(
    "python", "tgc_wsl_runner.py",
    package = "tidygenclust",
    mustWork = TRUE
  )
  .tgc_wsl_run_env(
    "ctidygenclust",
    "python",
    c(.tgc_wsl_path(runner), .tgc_wsl_path(config_path))
  )

  if (!file.exists(result_path)) {
    stop("The WSL runner did not produce a result manifest.", call. = FALSE)
  }
  manifest <- jsonlite::read_json(result_path, simplifyVector = TRUE)

  read_result_matrix <- function(entry) {
    .tgc_read_binary_matrix(
      file.path(job_dir, entry$file),
      entry$shape
    )
  }

  q_matrix <- read_result_matrix(manifest$Q)
  p_matrix <- if (!is.null(manifest$P)) {
    read_result_matrix(manifest$P)
  } else {
    NULL
  }
  cv_result <- if (!is.null(manifest$cv)) {
    list(avg = as.numeric(manifest$cv))
  } else {
    NULL
  }

  # Preserve all four return layouts: Q, Q/CV, Q/P, or Q/P/CV.
  no_freqs <- isTRUE(args$no_freqs)
  has_cv <- !is.null(args$cv)
  if (no_freqs && !has_cv) {
    return(q_matrix)
  }
  if (no_freqs && has_cv) {
    return(list(q_matrix, cv_result))
  }
  if (!has_cv) {
    return(list(q_matrix, p_matrix))
  }
  list(q_matrix, p_matrix, cv_result)
}

# Pass the public installer's source pins and choices to the bundled script.
# Environment creation, dependency locks and reuse are handled there.
.tgc_wsl_install_tools <- function(
  reset,
  fastmixture_hash,
  clumppling_hash,
  adamixture,
  adamixture_hash
) {
  if (!.tgc_wsl_available()) {
    stop(
      "Cannot start WSL distribution '", .tgc_wsl_distro(), "'. ",
      "Run 'wsl --list --verbose' in PowerShell and check its name. ",
      "Install WSL2 and Ubuntu if needed; see ",
      "vignette('windows', package = 'tidygenclust').",
      call. = FALSE
    )
  }

  installer <- system.file(
    "wsl", "install_tidygenclust.sh",
    package = "tidygenclust",
    mustWork = TRUE
  )
  status <- .tgc_wsl_call(c(
    "/bin/sh",
    .tgc_wsl_path(installer),
    .tgc_wsl_root(),
    fastmixture_hash,
    clumppling_hash,
    if (isTRUE(reset)) "true" else "false",
    if (isTRUE(adamixture)) "true" else "false",
    adamixture_hash
  ))
  if (!identical(as.integer(status), 0L)) {
    stop(
      "The tidygenclust WSL tool installation failed with exit status ",
      status, ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Query Linux package metadata; missing optional tools are reported as NA.
.tgc_wsl_package_version <- function(envname, package) {
  code <- paste0(
    "import importlib.metadata as m; ",
    "print(m.version(", encodeString(package, quote = "'"), "))"
  )
  out <- tryCatch(
    .tgc_wsl_run_env(
      envname, "python", c("-c", code),
      stdout = TRUE, stderr = TRUE
    ),
    error = function(e) NA_character_
  )
  if (length(out) == 0L || !is.null(attr(out, "status"))) {
    return(NA_character_)
  }
  trimws(out[[length(out)]])
}
