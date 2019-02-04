#' alpha_estimate
#'
#' Build objects of class alpha_estimate from their components. \code{alpha_estimate()} is a constructor method
#'
#' @param estimate The estimate
#' @param error The standard error in the estimate
#' @param estimand What is the estimate trying to estimate? (richness, Shannon...)
#' @param name The name of the method
#' @param interval An interval estimate
#' @param interval_type Type of interval estimate
#' @param type TODO(Amy): Deprecate?
#' @param model  What model is fit
#' @param warnings  Any warnings?
#' @param frequentist Logical. Frequentist or Bayesian?
#' @param parametric Logical. Parametric or not?
#' @param plot A ggplot associated with the estimate
#' @param reasonable Is the estimate likely to be reasonable?
#' @param other Any other relevant objects
#' @param ... Any other objects
#'
#' @return An object of class alpha_estimate
#'
#' @export
alpha_estimate <- function(estimate = NULL,
                           error = NULL,
                           estimand = NULL,
                           name = NULL,
                           interval = NULL,
                           interval_type = NULL,
                           type = NULL,
                           model = NULL,
                           warnings = NULL,
                           frequentist = NULL,
                           parametric = NULL,
                           plot = NULL,
                           reasonable = NULL,
                           other = NULL,
                           ...) {

  # if (is.null(ci) & !is.null(error)) {
  # # TODO need f0
  #   d <- exp(1.96*sqrt(log(1 + error^2 / f0)))
  #
  # }

  alpha_object <- list(estimate = estimate,
                       error = error,
                       estimand = tolower(estimand),
                       name = name,
                       interval = interval, #ifelse(is.na(estimate), c(NA, NA), interval),
                       interval_type = interval_type,
                       type = type,
                       model = model,
                       warnings = warnings,
                       frequentist = frequentist,
                       parametric = parametric,
                       plot = plot,
                       reasonable = reasonable,
                       other = other,
                       est = estimate,
                       seest = error,
                       ci = interval,
                       ...)


  class(alpha_object) <- append("alpha_estimate", class(alpha_object))

  return(alpha_object)
}

#' @export
print.alpha_estimate <- function(x, ...) {

  if (is.null(x$estimand)) {
    cat("Attempt to print estimate with unknown estimand")
  } else if (is.null(x$name)) {
    cat("Attempt to print estimate with unknown name")
  } else if (is.null(x$estimate)) {
    cat("Attempt to print estimate with unknown estimate")
  } else {
    cat(paste("Estimate of ", x$estimand,
              " from method ", x$name, ":\n", sep=""))
    cat(paste("  Estimate is ", round(x$estimate, ifelse(x$estimand == "richness", 0, 2)),
              " with standard error ", round(x$error, 2), "\n", sep=""))
    if (!is.null(x$interval)) {
      cat(paste("  Confidence interval: (",
                round(x$interval[1], ifelse(x$estimand == "richness", 0, 2)), ", ",
                round(x$interval[2], ifelse(x$estimand == "richness", 0, 2)), ")\n", sep=""))
    }
    if (!is.null(x$other$cutoff)) {
      cat(paste("  Cutoff: ", x$other$cutoff))
    }
    cat("\n")
  }
}

#' @export
summary.alpha_estimate <- function(object, ...) {
  # output just like a list

  # don't plot
  y <- object
  class(y) <- setdiff(class(y), "alpha_estimate")
  y$plot <- NULL
  print(y)
}


#' @export
plot.alpha_estimate <- function(x, ...) {
  x$plot
}




#' convert between different inputs for alpha-diversity estimates
#'
#' Inputs slated for development include phyloseq and otu_table
#'
#' @param input_data Supported types include filenames, data frames, matrices, vectors...
#'
#' @return Frequency count able
#'
#' @export
convert <- function(input_data) {

  if (class(input_data) == "character") {

    stop("breakaway no longer supports file paths as inputs")

  } else if ("data.frame" %in% class(input_data) |
             "matrix" %in% class(input_data)) {

    # determine if frequency count table or vector
    if (dim(input_data)[2] != 2) {
      stop("input_data is a data.frame or matrix but not a frequency count table.\n")
    }

    output_data <- input_data

  } else if (class(input_data) %in% c("numeric", "integer")) {

    # must be vector of counts
    if (isTRUE(all.equal(sum(input_data), 1)) &
        length(unique(input_data)) > 2) {
      stop("species richness estimates cannot accept relative abundances")
    }

    if (any(input_data %% 1 != 0)) {
      stop("input_data is a vector but not a vector of integers")
    }

    input_data_remove_zeros <- input_data[input_data != 0]
    output_data <- as.data.frame(table(input_data_remove_zeros))

  } else {
    stop(paste("Unsupported input type to function `convert`.",
               "You passed in an object of class", class(input_data)))
  }

  checked_output_data <- check_format(output_data)
  checked_output_data

}


#' Run some basic checks on a possible frequency count table
#'
#' @param output_data A matrix to test
#'
#' @return A checked frequency count table
check_format <- function(output_data) {

  if (length(output_data) <= 1) {
    warning("Input data to `check_format` is of length 1 or 0.")
    return(NULL)
  }

  if(!(class(output_data) %in% c("matrix", "data.frame"))) stop("Input should be a matrix or a data frame")

  if(length(dim(output_data)) != 2) stop("Input should have 2 columns")

  if(any(output_data[,2] %% 1 != 0)) stop("Second input column not integer-valued; should be counts")

  if(!all(rank(output_data[,1]) == 1:length(output_data[,1]))) warning("Frequency count format, right?")


  ## if table is used to create the frequency tables, the frequency index column is usually a factor, so fix this here
  if (class(output_data[,1]) == "factor") {
    as.integer(as.character(output_data[,1]))
  }

  # remove rows with zero
  output_data <- output_data[output_data[, 2] != 0, ]
  output_data <- data.frame(output_data)

  colnames(output_data) <- c("index", "frequency")

  if (!all(sort(output_data[, 1]) == output_data[, 1])) {
    stop("frequency counts not in order in `convert`")
  }

  output_data
}
