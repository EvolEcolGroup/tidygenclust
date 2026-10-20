# handles for the python modules
.py_rfastmixture <- NULL
.py_radamixture <- NULL
.onLoad <- function(...) {
  if (.Platform$OS.type == "windows") {
    # Windows uses separate WSL processes; do not initialise reticulate here.
    return(invisible(NULL))
  }

  if (reticulate::condaenv_exists("ctidygenclust")) {
    reticulate::use_condaenv("ctidygenclust", required = TRUE)
    .py_rfastmixture <<- reticulate::import_from_path(
      module = "py_rfastmixture",
      path = system.file("python", package = "tidygenclust"),
      delay_load = TRUE
    )
    .py_radamixture <<- reticulate::import_from_path(
      module = "py_radamixture",
      path = system.file("python", package = "tidygenclust"),
      delay_load = TRUE
    )
  } else {
    warning(
      "The conda environment 'ctidygenclust' does not exist. Install it ",
      "with 'tgc_tools_install()'"
    )
  }
}
.onAttach <- function(...) {
  if (.Platform$OS.type == "windows") {
    if (!nzchar(Sys.which("wsl.exe"))) {
      packageStartupMessage(
        "tidygenclust requires WSL2 on Windows. Install WSL, then run ",
        "'tgc_tools_install()'."
      )
    }
    return(invisible(NULL))
  }

  if (reticulate::condaenv_exists("ctidygenclust")) {
    reticulate::use_condaenv("ctidygenclust", required = TRUE)
  } else {
    packageStartupMessage(
      "The conda environment 'ctidygenclust' does not exist. Install it ",
      "with 'tgc_tools_install()'"
    )
  }
}
