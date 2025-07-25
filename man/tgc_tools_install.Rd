% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tgc_tools_install.R
\name{tgc_tools_install}
\alias{tgc_tools_install}
\title{Install tools for \code{tidygenclust}}
\usage{
tgc_tools_install(
  reset = FALSE,
  fastmixture_hash = "105eb99248d278cad320885190b919ad8a69be1b",
  clumppling_hash = "a4bf351037fb569e2c2cb83c603a1931606d4d40"
)
}
\arguments{
\item{reset}{a boolean used to reset the virtual environment. Only set
to TRUE if you have a broken virtual environment that you want to reset.}

\item{fastmixture_hash}{a string with the commit hash of the \code{fastmixture}
version to install. Default is the latest tested version.}

\item{clumppling_hash}{a string with the commit hash of the \code{clumppling}
version to install. Default is the latest tested version.}
}
\description{
\code{tidygenclust} relies on the \code{admxiture} and on python packages \code{fastmixture}
and \code{clumppling} for a number of
functionalities. We use \code{reticulate} to install them in conda
environments. As their dependencies are incompatible, we use two separate
conda environment, \code{ctidygenclust} (for \code{fastmixture} and \code{admixture}) and
\code{cclumppling} (for \code{clumppling}).
For each tool, default to the latest tested version of
these packages that have been tested to work with \code{tidegenclust}. It is
possible to provide a more recent github commit for a specific tool, but
this might lead to incompatibilities and errors.
}
