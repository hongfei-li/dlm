################################################################################
### a collection of utility functions
### version controlled by git
### author: Wenjie
### update date: 2018-09-30
################################################################################


### functions related with I/O of files

##' Input and Output Direcotory
##'
##' @name io_dir
##'
##' @param path A character vector of length one for specifying working
##'     directory.
##' @param verbose A logical vector of length one.  Suppress verbose messages
##'     if \code{FALSE}.
##' @param ... Other arguments for future usage.
NULL

##' @rdname io_dir
##' @export
input_dir <- function(path = ".", ...)
{
    if (! dir.exists(path))
        stop("The input directory does not exist!", call. = FALSE)
    path
}

##' @rdname io_dir
##' @export
output_dir <- function(path = ".", verbose = TRUE, ...)
{
    if (! dir.exists(path)) {
        dir.create(path, recursive = TRUE)
        if (verbose) {
            message(sprintf("created %s", path))
        }
    }
    path
}


##' Check, Install and Attach Multiple R packages Specified
##'
##' The function is a simple wrapper that aims to simplify the process of
##' attaching multiple packages that are not necessarily available.  It first
##' checks whether the packages given were already installed and would try to
##' install them from specified repository if needed.  At last, it attachs all
##' the needed packages to the search path.
##'
##' See the source code for details.
##'
##' @param pkgs A character vector specifying the packages needed to reproduce
##'     this document.
##' @param repos A character vector for the base URL(s) of the repositories
##'     containing the source.  The default mirror is
##'     \code{"https://cloud.r-project.org"} and will be used if \code{repos}
##'     is not set.
##' @param ... Other arguments passed to function
##'     \code{\link[utils]{install.packages}}.
##'
##' @return \code{NULL} invisibly.
##' @examples
##' \dontrun{
##' need.pacakges(c("data.table", "ggplot2"))
##' }
##' @importFrom utils installed.packages install.packages
##' @export
need.packages <- function(pkgs, repos = getOption("repos"), ...)
{
    new.pkgs <- pkgs[! (pkgs %in% installed.packages()[, "Package"])]
    if (length(new.pkgs) > 0) {
        if (is.null(repos) || repos == "@CRAN@") {
            repos <- "https://cloud.r-project.org"
        }
        install.packages(pkgs = new.pkgs, repos = repos, ...)
    }
    sapply(pkgs, function(a) {
        suppressMessages(require(a, character.only = TRUE))
    })
    invisible()
}
