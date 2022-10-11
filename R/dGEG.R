#' The Generalised exponential-Gaussian distribution
#'
#' @description
#' Density, distribution function, quantile function,
#' random generation and hazard function for the Generalised exponential-Gaussian distribution with
#' parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}.
#'
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param tau parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#'
#' @details
#' The Generalised exponential-Gaussian with parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}
#' has density given by
#'
#' \eqn{f(x | \mu, \sigma, \nu, \tau) = \frac{\tau}{\nu} \exp(w) \Phi \left( z - \frac{\sigma}{\nu} \right) \left[ \Phi(z) - \exp(w)  \Phi \left( z - \frac{\sigma}{\nu} \right) \right]^{\tau-1}}
#'
#' for \eqn{-\infty < x < \infty}. With \eqn{w=\frac{\mu-x}{\nu} + \frac{\sigma^2}{2\nu^2}} and \eqn{z=\frac{x-\mu}{\sigma}}.
#'
#' @return
#' \code{dGEG} gives the density, \code{pGEG} gives the distribution
#' function, \code{qGEG} gives the quantile function, \code{rGEG}
#' generates random deviates.
#'
#' @example examples/examples_dGEG.R
#'
#' @importFrom stats pnorm
#' @export
dGEG <- function(x, mu=0, sigma=1, nu=1, tau=1, log=FALSE){
  if (any(sigma <= 0))
    stop("parameter sigma has to be positive!")
  if (any(nu <= 0))
    stop("parameter nu has to be positive!")
  if (any(tau <= 0))
    stop("parameter tau has to be positive!")
  z <- (x-mu)/sigma
  w <- (mu-x)/nu + sigma^2/(2*nu^2)
  res <- log(tau) - log(nu) + w
  res <- res + log(pnorm(z - sigma/nu))
  res <- res + (tau-1) * log(pnorm(z) - exp(w) * pnorm(z - sigma/nu))
  if(log)
    return(res)
  else
    return(exp(res))
}
#' @importFrom stats pnorm
#' @export
#' @rdname dGEG
pGEG <- function(q, mu=0, sigma=1, nu=1, tau=1, lower.tail=TRUE, log.p=FALSE){
  if (any(sigma <= 0))
    stop("parameter sigma has to be positive!")
  if (any(nu <= 0))
    stop("parameter nu has to be positive!")
  if (any(tau <= 0))
    stop("parameter tau has to be positive!")
  z <- (q-mu)/sigma
  w <- (mu-q)/nu + sigma^2/(2*nu^2)
  cdf <- tau * log(pnorm(z) - exp(w) * pnorm(z - sigma/nu))
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- cdf
  else cdf <- exp(cdf)
  cdf
}
#' @importFrom gamlss.dist qexGAUS
#' @export
#' @rdname dGEG
qGEG <- function(p, mu=0, sigma=1, nu=1, tau=1){
  if (any(sigma <= 0))
    stop("parameter sigma has to be positive!")
  if (any(nu <= 0))
    stop("parameter nu has to be positive!")
  if (any(tau <= 0))
    stop("parameter tau has to be positive!")
  #library(gamlss)
  return(qexGAUS(p=p^(1/tau), mu=mu, sigma=sigma, nu=nu))
}

#' @importFrom stats runif
#' @importFrom gamlss.dist qexGAUS
#' @export
#' @rdname dGEG
rGEG <- function(n=1, mu=0, sigma=1, nu=1, tau=1){
  if (any(sigma <= 0))
    stop("parameter sigma has to be positive!")
  if (any(nu <= 0))
    stop("parameter nu has to be positive!")
  if (any(tau <= 0))
    stop("parameter tau has to be positive!")
  u <- runif(n, 0, 1)^(1/tau)
  #library(gamlss)
  x <- qexGAUS(p=u, mu=mu, sigma=sigma, nu=nu)
  return(x)
}
