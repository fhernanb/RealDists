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
