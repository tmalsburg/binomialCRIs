
library(ggplot2)

#' Plots posterior and HPD interval for the probability parameter
#' underlying a binomially distributed outcome.
#'
#' @param n.successes number of successes.
#' @param n.trials total number of trials.
#' @param prior.alpha alpha parameter of the Beta distribution
#'   defining the prior.  The default values alpha=1 and beta=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior.beta beta parameter of the Beta distribution defining
#'   the prior.  The default values alpha=1 and beta=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @param prob the size of the HPDI interval with 0 < prob <= 1.
#' @return ggplot2 object.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
plot.binomial.hpdi <- function(n.successes, n.trials, prior.alpha=1, prior.beta=1, prob=0.9) {
  
  interval <- binomial.hpdi(n.successes, n.trials, prior.alpha, prior.beta, prob)
  
  plot.binomial.cri(interval, n.successes, n.trials, prior.alpha, prior.beta) +
    labs(caption=paste("Data: successes=", n.successes, ", trials=", n.trials, "\nPrior: Beta distribution with α=", prior.alpha, ", β=", prior.beta, "\n", 100*prob, "% HPD interval: ", round(interval[1], 2), "--", round(interval[2], 2), sep=""))
  
}

#' Plots posterior and percentile interval for the probability
#' parameter underlying a binomially distributed outcome.
#'
#' @param n.successes number of successes.
#' @param n.trials total number of trials.
#' @param prior.alpha alpha parameter of the Beta distribution
#'   defining the prior.  The default values alpha=1 and beta=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior.beta beta parameter of the Beta distribution defining
#'   the prior.  The default values alpha=1 and beta=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @param prob the size of the percentile interval with 0 <= prob <= 1.
#' @return ggplot2 object.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
plot.binomial.pi <- function(n.successes, n.trials, prior.alpha=1, prior.beta=1, prob=0.9) {
  
  interval <- binomial.pi(n.successes, n.trials, prior.alpha, prior.beta, prob)
  
  plot.binomial.cri(interval, n.successes, n.trials, prior.alpha, prior.beta) +
    labs(caption=paste("Data: successes=", n.successes, ", trials=", n.trials, "\nPrior: Beta distribution with α=", prior.alpha, ", β=", prior.beta, "\n", 100*prob, "% precentile interval: ", round(interval[1], 2), "--", round(interval[2], 2), sep=""))
  
}

# Plots posterior distribution (Beta) with given interval shaded:
plot.binomial.cri <- function(interval, n.successes, n.trials, prior.alpha=1, prior.beta=1) {
  
  alpha <- n.successes + prior.alpha
  beta  <- n.trials - n.successes + prior.beta

  p_grid <- seq(0, 1, length.out=1000)
  dens   <- dbeta(p_grid, alpha, beta)
  
  # For plotting the interval:
  x.interval <- p_grid[p_grid >= interval[1] & p_grid <= interval[2]]
  y.interval <-   dens[p_grid >= interval[1] & p_grid <= interval[2]]
  x.interval <- x.interval[c(1, 1:length(x.interval), length(x.interval))]
  y.interval <- c(0, y.interval, 0)
  
  ggplot(mapping=aes(x=p_grid, y=dens)) +
    geom_line() +
    geom_polygon(aes(x = x.interval, y = y.interval), alpha=0.5) +
    labs(x="θ",
         y="Posterior density")
  
}

#' Calculates the HPD interval for the probability parameter
#' underlying a binomially distributed outcome.
#'
#' @param n.successes number of successes.
#' @param n.trials total number of trials.
#' @param prior.alpha alpha parameter of the Beta distribution
#'   defining the prior.  The default values alpha=1 and beta=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior.beta beta parameter of the Beta distribution defining
#'   the prior.  The default values alpha=1 and beta=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @param prob the size of the HPDI interval with 0 < prob <= 1.
#' @return vector containing the endpoints of the HPDI.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
binomial.hpdi <- function(n.successes, n.trials, prior.alpha=1, prior.beta=1, prob=0.9) {

  stopifnot(n.successes == as.integer(n.successes))
  stopifnot(n.trials == as.integer(n.trials))
  stopifnot(n.trials > 0)
  stopifnot(n.successes > 0)
  stopifnot(n.successes <= n.trials)
  stopifnot(prior.alpha > 0)
  stopifnot(prior.beta > 0)
  stopifnot(prob > 0 & prob <= 1)
  
  if(prob==1) {
    x <- c(0, 1)
  } else {
    # Parameters for posterior distribution:
    alpha <- n.successes + prior.alpha
    beta  <- n.trials - n.successes + prior.beta
    
    f <- function(l) {
      qbeta(l+prob, alpha, beta) - qbeta(l, alpha, beta)
    }
    
    x <- optimize(f, c(0, 1-prob), maximum=FALSE)$minimum
    x <- c(qbeta(x, alpha, beta), qbeta(x + prob, alpha, beta))
  }
  
  names(x) <- c(paste("|", prob, sep=""), paste(prob, "|", sep=""))
  x
}

#' Calculates the percentile interval for the probability parameter
#' underlying a binomially distributed outcome.
#'
#' @param n.successes number of successes.
#' @param n.trials total number of trials.
#' @param prior.alpha alpha parameter of the Beta distribution
#'   defining the prior.  The default values alpha=1 and beta=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior.beta beta parameter of the Beta distribution defining
#'   the prior.  The default values alpha=1 and beta=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @param prob the size of the percentile interval with 0 <= prob <= 1.
#' @return vector containing the endpoints of the percentile interval.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
binomial.pi <- function(n.successes, n.trials, prior.alpha=1, prior.beta=1, prob=0.9) {

  stopifnot(n.successes == as.integer(n.successes))
  stopifnot(n.trials == as.integer(n.trials))
  stopifnot(n.trials > 0)
  stopifnot(n.successes > 0)
  stopifnot(n.successes <= n.trials)
  stopifnot(prior.alpha > 0)
  stopifnot(prior.beta > 0)
  stopifnot(prob >= 0 & prob <= 1)

  # Parameters for posterior distribution:
  alpha <- n.successes + prior.alpha
  beta  <- n.trials - n.successes + prior.beta
  
  x <- c(qbeta((1-prob)/2,        alpha, beta),
         qbeta((1-prob)/2 + prob, alpha, beta))
  
  names(x) <- c(paste("|", prob, sep=""), paste(prob, "|", sep=""))
  x
}

#' Calculates the posterior probability that the probability parameter
#' underlying a binomially distributed outcome is larger than some
#' threshold.
#'
#' @param n.successes number of successes.
#' @param n.trials total number of trials.
#' @param prior.alpha alpha parameter of the Beta distribution
#'   defining the prior.  The default values alpha=1 and beta=1 define
#'   a flat prior assigning equal probability density to all possible
#'   parameter values.
#' @param prior.beta beta parameter of the Beta distribution defining
#'   the prior.  The default values alpha=1 and beta=1 define a flat
#'   prior assigning equal probability density to all possible
#'   parameter values.
#' @param prob the threshold.
#' @return Posterior probability that the parameter value is larger
#'   than the threshold \code{prob}.
#' @author Titus von der Malsburg <malsburg@uni-potsdam.de>
binomial.pgt <- function(n.successes, n.trials, prior.alpha=1, prior.beta=1, prob=0.5) {
  
  stopifnot(n.successes == as.integer(n.successes))
  stopifnot(n.trials == as.integer(n.trials))
  stopifnot(n.trials > 0)
  stopifnot(n.successes > 0)
  stopifnot(n.successes <= n.trials)
  stopifnot(prior.alpha > 0)
  stopifnot(prior.beta > 0)
  stopifnot(prob >= 0 & prob <= 1)
  
  # Parameters for posterior distribution:
  alpha <- n.successes + prior.alpha
  beta  <- n.trials - n.successes + prior.beta

  1 - pbeta(prob, alpha, beta)
  
}

