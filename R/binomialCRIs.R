#' Given data and a prior, the functions in the package calculate and
#' plot credible intervals for the probability of success.
#'
#' @name binomialCRIs
#' @docType package
#' @title Calculate and Plot Binomial Credible Intervals
#' @author Titus von der Malsburg \email{malsburg@@uni-potsdam.de}
#' @keywords Bayesian inference, binomial distribution, credible intervals
#' @seealso \code{\link{binomial_hpdi}}, \code{\link{binomial_pi}}, \code{\link{binomial_prob}}, \code{\link{plot_binomial_hpdi}}, \code{\link{plot_binomial_pi}}, \code{\link{plot_binomial_cri}}

NULL

check.parameters <- function(n_successes, n_trials, prior_shape1, prior_shape2) {
  stopifnot(n_successes == as.integer(n_successes))
  stopifnot(n_trials == as.integer(n_trials))
  stopifnot(0 <= n_successes & n_successes <= n_trials)
  stopifnot(0 < n_trials)
  stopifnot(prior_shape1 > 0)
  stopifnot(prior_shape2 > 0)
}

#' Calculates the HPD interval for the probability parameter
#' underlying a binomially distributed outcome.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob The size of the HPDI interval with 0 < prob <= 1.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return A vector containing the endpoints of the HPDI.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' # 6/9 successes with flat prior:
#' binomial_hpdi(6, 9, 0.8)
#' 
#' # 6/9 successes with prior assuming one prior success and one failure:
#' binomial_hpdi(6, 9, 0.8, 2, 2)
#' 
#' # Equivalently:
#' binomial_hpdi(7, 11, 0.8)
binomial_hpdi <- function(n_successes, n_trials, prob, prior_shape1=1, prior_shape2=1) {

  check.parameters(n_successes, n_trials, prior_shape1, prior_shape2)
  stopifnot(0 < prob & prob <= 1)
  
  if(prob==1) {
    x <- c(0, 1)
  } else {
    # Parameters for posterior distribution:
    shape1 <- n_successes + prior_shape1
    shape2 <- n_trials - n_successes + prior_shape2
    
    f <- function(l) {
      stats::qbeta(l+prob, shape1, shape2) - stats::qbeta(l, shape1, shape2)
    }
    
    x <- stats::optimize(f, c(0, 1-prob), maximum=FALSE)$minimum
    x <- c(stats::qbeta(x, shape1, shape2), stats::qbeta(x + prob, shape1, shape2))
  }
  
  names(x) <- c(paste("|", prob, sep=""), paste(prob, "|", sep=""))
  x
}

#' Calculates the percentile interval for the probability parameter
#' underlying a binomially distributed outcome.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob The size of the percentile interval with 0 <= prob <= 1.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return A vector containing the endpoints of the percentile interval.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' # 6/9 successes with flat prior:
#' binomial_pi(6, 9, 0.8)
#' 
#' # 6/9 successes with prior assuming one prior success and one failure:
#' binomial_pi(6, 9, 0.8, 2, 2)
#' 
#' # Equivalently:
#' binomial_pi(7, 11, 0.8)
binomial_pi <- function(n_successes, n_trials, prob, prior_shape1=1, prior_shape2=1) {

  check.parameters(n_successes, n_trials, prior_shape1, prior_shape2)
  stopifnot(0 <= prob & prob <= 1)

  # Parameters for posterior distribution:
  shape1 <- n_successes + prior_shape1
  shape2 <- n_trials - n_successes + prior_shape2
  
  x <- c(stats::qbeta((1-prob)/2,        shape1, shape2),
         stats::qbeta((1-prob)/2 + prob, shape1, shape2))
  
  names(x) <- c(paste("|", prob, sep=""), paste(prob, "|", sep=""))
  x
}

#' Calculates the posterior probability that the probability parameter
#' underlying a binomially distributed outcome is in a specified
#' interval.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob_lower The lower end point of the interval.
#' @param prob_upper The upper end point of the interval.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return The posterior probability that the parameter value lies in the
#'   specified interval.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' # Probability of parameter being larger than 0.5 after seeing 6/9
#' # successes with flat prior:
#' binomial_prob(6, 9, 0.5)
#' 
#' # Probability of parameter being smaller than 0.5 after seeing 6/9
#' # successes with prior assuming one earlier success and one
#' # failure:
#' binomial_prob(6, 9, prob_upper=0.5, prior_shape1=2, prior_shape2=2)
#'
#' # Probability of parameter being larger than 0.5 and smaller than
#' # 0.75 after seeing 6/9 successes with flat prior:
#' binomial_prob(6, 9, 0.5, 0.75)
binomial_prob <- function(n_successes, n_trials, prob_lower=0, prob_upper=1, prior_shape1=1, prior_shape2=1) {
  
  check.parameters(n_successes, n_trials, prior_shape1, prior_shape2)
  stopifnot(0 <= prob_lower & prob_lower <= prob_upper & prob_upper <= 1)
  
  # Parameters for posterior distribution:
  shape1 <- n_successes + prior_shape1
  shape2 <- n_trials - n_successes + prior_shape2

  stats::pbeta(prob_upper, shape1, shape2) - stats::pbeta(prob_lower, shape1, shape2)
  
}

#' Plots the posterior distribution and specified HPD interval for the
#' probability parameter underlying a binomially distributed outcome.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob The size of the HPDI interval with 0 < prob <= 1.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return A ggplot2 object.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' plot_binomial_hpdi(6, 9, 0.5)
plot_binomial_hpdi <- function(n_successes, n_trials, prob, prior_shape1=1, prior_shape2=1) {
  
  interval <- binomial_hpdi(n_successes, n_trials, prob, prior_shape1, prior_shape2)
  
  plot_binomial_cri(n_successes, n_trials, interval[1], interval[2], prior_shape1, prior_shape2) +
    ggplot2::labs(caption=paste("Data: successes=", n_successes, ", trials=",
                              n_trials, "\nPrior: Beta distribution with \u3b1=", prior_shape1,
                                ", \u3b2=", prior_shape2, "\n", 100*prob, "% HPD interval: ",
                                round(interval[1], 2), "\u2013", round(interval[2], 2), sep=""))
}

#' Plots the posterior distribution and specified percentile interval
#' for the probability parameter underlying a binomially distributed
#' outcome.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob The size of the percentile interval with 0 <= prob <= 1.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return A ggplot2 object.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' plot_binomial_pi(6, 9, 0.5)
plot_binomial_pi <- function(n_successes, n_trials, prob, prior_shape1=1, prior_shape2=1) {
  
  interval <- binomial_pi(n_successes, n_trials, prob, prior_shape1, prior_shape2)
  
  plot_binomial_cri(n_successes, n_trials, interval[1], interval[2], prior_shape1, prior_shape2) +
    ggplot2::labs(caption=paste("Data: successes=", n_successes, ", trials=",
                                n_trials, "\nPrior: Beta distribution with \u3b1=", prior_shape1,
                                ", \u3b2=", prior_shape2, "\n", 100*prob, "% precentile interval: ",
                                round(interval[1], 2), "\u2013", round(interval[2], 2), sep=""))
}

#' Plots the posterior distribution and a specified interval (if
#' defined) for the probability parameter underlying a binomially
#' distributed outcome.
#'
#' @param n_successes The number of successes.
#' @param n_trials The total number of trials.
#' @param prob_lower The lower end point of the interval. Defaults to
#'   0 if \code{prob_upper} is non-null.
#' @param prob_upper The upper end point of the interval. Defaults to
#'   1 if \code{prob_lower} is non-null.
#' @param prior_shape1 The shape1 parameter of the Beta distribution
#'   defining the prior.  The default values shape1=1 and shape2=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior_shape2 The shape2 parameter of the Beta distribution defining
#'   the prior.  The default values shape1=1 and shape2=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @return A ggplot2 object.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
#' @export
#' @examples
#' plot_binomial_cri(6, 9, 0.5, 0.75)
plot_binomial_cri <- function(n_successes, n_trials, prob_lower=NULL, prob_upper=NULL, prior_shape1=1, prior_shape2=1) {
  
  check.parameters(n_successes, n_trials, prior_shape1, prior_shape2)
  stopifnot(0 <= prob_lower & prob_lower <= prob_upper & prob_upper <= 1)

  if (!is.null(prob_lower) & is.null(prob_upper)) prob_upper <- 1
  if (!is.null(prob_upper) & is.null(prob_lower)) prob_lower <- 0
  
  shape1 <- n_successes + prior_shape1
  shape2 <- n_trials - n_successes + prior_shape2

  p_grid <- seq(0, 1, length.out=1000)
  dens   <- stats::dbeta(p_grid, shape1, shape2)
  
  # For plotting the interval:
  x_interval <- p_grid[p_grid >= prob_lower & p_grid <= prob_upper]
  y_interval <-   dens[p_grid >= prob_lower & p_grid <= prob_upper]
  x_interval <- x_interval[c(1, 1:length(x_interval), length(x_interval))]
  y_interval <- c(0, y_interval, 0)

  if (!is.null(prob_lower)) {
    p <- binomial_prob(n_successes, n_trials, prob_lower, prob_upper, prior_shape1, prior_shape2)
    p <- paste("\nPr(\u3b8 in [", prob_lower, ",", prob_upper, "]): ", format(p, digits=2), sep="")
  } else { 
    p <- ""
  }
    
  ggplot2::ggplot(mapping=ggplot2::aes(x=p_grid, y=dens)) +
    ggplot2::geom_line() +
    ggplot2::geom_polygon(ggplot2::aes(x=x_interval, y=y_interval), alpha=0.5) +
    ggplot2::labs(x="\u3b8",
                  y="Posterior density",
                  caption=paste("Data: successes=", n_successes, ", trials=",
                                n_trials, "\nPrior: Beta distribution with \u3b1=", prior_shape1,
                                ", \u3b2=", prior_shape2, p, sep=""))
  
}

