
# Example 1 - without covariates ------------------------------------------

n <- 100

# The true parameters are:
true_mu <- -5
true_si <- 4
true_nu <- 2.5
true_ta <- 1
true_theta <- c(true_mu, true_si, true_nu, true_ta)

# Graphing the pdf
curve(dGEG(x, mu=true_mu, sigma=true_si, nu=true_nu, tau=true_ta),
      ylab="Density", xlab="X", las=1,
      from=-25, to=25, lwd=3, col="tomato")

# Simulating a random sample
y <- rGEG(n=n, mu=true_mu, sigma=true_si, nu=true_nu, tau=true_ta)

# Estimating paramaters
library(gamlss)
mod <- gamlss(y ~ 1, family=GEG,
              control=gamlss.control(n.cyc=1000, trace=TRUE))

# Vector with the estimated results
res <- c(mu_hat=coef(mod, what="mu"),
         sigma_hat=exp(coef(mod, what="sigma")),
         nu_hat=exp(coef(mod, what="nu")),
         tau_hat=exp(coef(mod, what="tau")))

# Comparing true vector and estimated vector
round(cbind(true_theta, with_GEG=res), digits=2)

# Histogram, estimated density and true density
truehist(y, ylab="Density", col="gray", las=1)

curve(dGEG(x, mu=res[1], sigma=res[2], nu=res[3], tau=res[4]),
      add=TRUE, col="blue2", lwd=2)

curve(dGEG(x, mu=true_theta[1], sigma=true_theta[2],
           nu=true_theta[3], tau=true_theta[4]),
      add=TRUE, col="green4", lwd=2)

legend("topright", lwd=2, bty="n",
       legend=c("with GEG", "true density"),
       col=c("blue2", "green4"))


# Example 2 - with covariates ---------------------------------------------
n <- 500

# The true parameters are:
b0_mu <- -1
b1_mu <-  2

b0_sigma <- -2
b1_sigma <-  4

true_nu <- 0.5
true_ta <- 0.75

# The true theta vector
true_theta <- c(b0_mu, b1_mu, b0_sigma, b1_sigma, true_nu, true_ta)

# Simulating covariates
x1 <- runif(n, min=0.49, max=0.51)
x2 <- runif(n, min=0.49, max=0.51)

# Simulating a random sample
y <- rGEG(n=n,
          mu    =     b0_mu    + b1_mu    * x1,
          sigma = exp(b0_sigma + b1_sigma * x2),
          nu    = true_nu,
          tau   = true_ta)

# The dataframe
datos <- data.frame(y=y, x1=x1, x2=x2)

# Estimating paramaters
# Using gamlss with our proposal
mod <- gamlss(y ~ x1,
              sigma.fo = ~ x2,
              family=GEG,
              control=gamlss.control(n.cyc=10000, trace=TRUE))

# To obtain the estimated parameters
param <- unlist(coefAll(mod))

# Comparing true vector and estimated vector
res <- cbind(true_theta, with_gamlss=c(param[1:4], exp(param[5:6])))
round(res, digits=2)

