

#' The transformed weighted linear regression estimator for species richness estimation
#'
#' This function implements the transformed version of the species richness
#' estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).
#'
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#' @param cutoff Maximum frequency count to use
#' @param print Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#'
#' @note While robust to many different structures, model is almost always
#' misspecified. The result is usually implausible diversity estimates with
#' artificially low standard errors. Extreme caution is advised.
#' @author Amy Willis
#' @seealso \code{\link{breakaway}}; \code{\link{apples}};
#' \code{\link{wlrm_untransformed}}
#' @references Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size
#' estimation based upon ratios of recapture probabilities. \emph{Annals of
#' Applied Statistics}, \bold{5}.
#' @keywords diversity models
#'
#' @import ggplot2
#' @import stats
#'
#' @examples
#'
#' wlrm_transformed(apples)
#' wlrm_transformed(apples, plot = FALSE, print = FALSE, answers = TRUE)
#'
#' @export

wlrm_transformed <- function(input_data, ...) {
  breakaway::wlrm_transformed(input_data, ...)
}


#' The untransformed weighted linear regression estimator for species richness estimation
#'
#' This function implements the untransformed version of the species richness
#' estimation procedure outlined in Rocchetti, Bunge and Bohning (2011).
#'
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#' @param cutoff Maximum frequency count to use
#' @param print Deprecated; only for backwards compatibility
#' @param plot Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @note This estimator is based on the negative binomial model and for that
#' reason generally produces poor fits to microbial data. The result is usually
#' artificially low standard errors. Caution is advised.
#'
#' @author Amy Willis
#' @seealso \code{\link{breakaway}}; \code{\link{apples}};
#' \code{\link{wlrm_transformed}}
#' @references Rocchetti, I., Bunge, J. and Bohning, D. (2011). Population size
#' estimation based upon ratios of recapture probabilities. \emph{Annals of
#' Applied Statistics}, \bold{5}.
#' @keywords diversity models
#' @importFrom graphics legend points
#' @importFrom utils head read.table
#'
#' @import ggplot2
#' @import stats
#'
#' @examples
#'
#' wlrm_untransformed(apples)
#'
#' @export

wlrm_untransformed <- function(input_data, ...) {
  breakaway::wlrm_untransformed(input_data, ...)
}


#' Chao1 species richness estimator
#'
#' This function implements the Chao1 richness estimate, which is often
#' mistakenly referred to as an index.
#'
#'
#' @param input_data An input type that can be processed by \code{convert()} or a \code{phyloseq} object
#' @param output Deprecated; only for backwards compatibility
#' @param answers Deprecated; only for backwards compatibility
#'
#' @return An object of class \code{alpha_estimate}, or \code{alpha_estimates} for \code{phyloseq} objects
#' @note The authors of this package strongly discourage the use of this
#' estimator.  It is only valid when you wish to assume that every taxa has
#' equal probability of being observed. You don't really think that's possible,
#' do you?
#' @author Amy Willis
#' @examples
#'
#'
#' chao1(apples)
#'
#'
#' @export
#'

chao1 <- function(input_data, ...) {
  breakaway::chao1(input_data, ...)
}


