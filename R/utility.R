##########################################################################
## Utility Functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


# takes a symmetric matrix x and returns lower diagonal
# note: does not check for symmetry
#
# ADM 4/18/2003

#' Extract Lower Triangular Elements from a Symmetric Matrix
#'
#' This function takes a symmetric matrix and extracts a list of all lower
#' triangular elements.
#'
#' This function checks to make sure the matrix is square, but it does not
#' check for symmetry (it just pulls the lower triangular elements).  The
#' elements are stored in column major order.  The original matrix can be
#' restored using the \code{xpnd} command.
#'
#' @param x A symmetric matrix.
#'
#' @return A list of the lower triangular elements.
#'
#' @export
#'
#' @seealso \code{\link{xpnd}}
#'
#' @keywords manip
#'
#' @examples
#'
#'    symmat <- matrix(c(1,2,3,4,2,4,5,6,3,5,7,8,4,6,8,9),4,4)
#'    vech(symmat)
#'
"vech" <-
  function (x) {
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
      stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
  }

# takes vector x and returns an nrow times nrow symmetric matrix
# this will recycle the elements of x as needed to fill the matrix
#
# ADM 4/18/2003
# ADM 11/13/2003 [bug fix]
# ADM 1/25/2006 [patch to automatically compute matrix size]

#' Expand a Vector into a Symmetric Matrix
#'
#' This function takes a vector of appropriate length (typically created using
#' \code{vech}) and creates a symmetric matrix.
#'
#' This function is particularly useful when dealing with variance covariance
#' matrices. Note that R stores matrices in column major order, and that the
#' items in \code{x} will be recycled to fill the matrix if need be.
#'
#' The number of rows can be specified or automatically computed from the
#' number of elements in a given object via \eqn{(-1 + \sqrt{(1 + 8 *
#' length(x))}) / 2}.
#'
#' @param x A list of elements to expand into symmetric matrix.
#'
#' @param nrow The number of rows (and columns) in the returned matrix.  Look
#' into the details.
#'
#' @export
#'
#' @return An \eqn{(nrows \times nrows)} symmetric matrix.
#'
#' @seealso \code{\link{vech}}
#'
#' @keywords manip
#'
#' @examples
#'
#'   xpnd(c(1,2,3,4,4,5,6,7,8,9),4)
#'   xpnd(c(1,2,3,4,4,5,6,7,8,9))
#'
"xpnd" <-
  function (x, nrow = NULL) {
    dim(x) <- NULL
    if(is.null(nrow)) nrow <- (-1 + sqrt(1 + 8 * length(x))) / 2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag=TRUE)] <- 0
    output <- output + t(hold)
    return(output)
  }
