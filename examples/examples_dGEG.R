## The probability density function
curve(dGEG(x, mu=7, sigma=0.25, nu=0.05, tau=1),
      from=0, to=10, col="red", las=1, ylab="f(x)")

curve(dGEG(x, mu=7, sigma=0.25, nu=0.05, tau=0.5),
      add=TRUE, col="blue3")

curve(dGEG(x, mu=7, sigma=0.25, nu=0.05, tau=0.1),
      add=TRUE, col="green4")

legend("topleft", col=c("red", "blue3", "green4"), lty=1, bty="n",
       legend=c("mu=7, sigma=0.25, nu=0.05, tau=1",
                "mu=7, sigma=0.25, nu=0.05, tau=0.5",
                "mu=7, sigma=0.25, nu=0.05, tau=0.1"))

## The cumulative distribution function
curve(pGEG(x, mu=-5, sigma=0.25, nu=0.05, tau=1), from=-10, to=10,
      col="red", las=1, ylab="F(x)")
curve(pGEG(x, mu=-2, sigma=0.25, nu=0.05, tau=0.5),
      add=TRUE, col="blue3", las=1)

## The quantile function
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qGEG(p, mu=7, sigma=0.25, nu=0.05, tau=1), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pGEG(x, mu=7, sigma=0.25, nu=0.05, tau=1), from=0, add=TRUE, col="red")

## The random function
hist(rGEG(n=1000, mu=7, sigma=0.25, nu=0.05, tau=1), freq=FALSE, xlab="x",
     ylim=c(0, 2), las=1, main="")
curve(dGEG(x, mu=7, sigma=0.25, nu=0.05, tau=1), add=TRUE, col="red")
