# handles for the python modules
.py_rclumppling <- NULL
.py_rfastmixture <- NULL
# .onLoad <- function(...) {
#   reticulate::use_virtualenv("rclumppling", required = FALSE)
#   .py_rclumppling <<- reticulate::import_from_path(module = "py_rclumppling",
#                                                   path = system.file("python",
#                                                                      package = "rclumppling"),
#                                                   delay_load = TRUE)
# }
.onLoad <- function(...) {
  .py_rclumppling <<- reticulate::import_from_path(
    module = "py_rclumppling",
    path = system.file("python", package = "tidygenclust"),
    delay_load = TRUE
  )
  .py_rfastmixture <<- reticulate::import_from_path(
    module = "py_rfastmixture",
    path = system.file("python", package = "tidygenclust"),
    delay_load = TRUE
  )
}



# .onLoad <- function(...) {
#   .py_rfastmixture <<- reticulate::import_from_path(
#     module = "py_rfastmixture",
#     path = system.file("python", package = "tidygenclust"),
#     delay_load = TRUE,
#     convert = FALSE # do not convert python objects to R objects automatically
#   )
# }

.onAttach <- function (...) {
  if (reticulate::condaenv_exists("ctidygenclust")) {
    reticulate::use_condaenv("ctidygenclust", required = FALSE)
  } else {
    packageStartupMessage(
      "The conda environment 'ctidygenclust' does not exist. Install it with tgc_tools_install()"
    )
  }
}
