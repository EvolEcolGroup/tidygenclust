# handles for the python modules
.py_rfastmixture <- NULL
.onLoad <- function(...) {
  .py_rfastmixture <<- reticulate::import_from_path(
    module = "py_rfastmixture",
    path = system.file("python", package = "tidygenclust"),
    delay_load = TRUE
  )
}
.onAttach <- function(...) {
  if (reticulate::condaenv_exists("ctidygenclust")) {
    reticulate::use_condaenv("ctidygenclust", required = FALSE)
  } else {
    packageStartupMessage(
      "The conda environment 'ctidygenclust' does not exist. Install it ",
      "with 'tgc_tools_install()'"
    )
  }
}
