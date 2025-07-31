preds <- rast("output/k13/gneiting_sparse/preds_all.grd")
ras <- rast("output/k13/surveillance_effort_k13.grd")

pwr_binom <- function(n, p){
  # probability k > 0 given n, p
  1 - (1 - p)^n
}

p <- preds$`2022_post_median`

# from surveillance_effort
gf <- focalMat(ras, 1, "Gauss")
n <- aggregate(ras, fact = 4, sum, na.rm = TRUE) %>%
  focal(gf, na.rm=TRUE) %>%
  project(p) %>%
  mask(p)

plot(n)

# aggregate makes more sense here ....
# n <- test_dens %>%
#   aggregate(fact = 4, sum, na.rm=TRUE) %>%
#   project(p) %>%
#   mask(p)
# n <- n + 1

# could also use extant test_dens + 1?

# shouldn't be a problem ...?
n[n < 1] <- NA

pwr <- pwr_binom(n, p)

plot(p, main="Estimated prevalence")
plot(sqrt(n), main="Surv intensity")
plot(pwr, main="Power:\n 'Can we detect if prevalence > 0 | prevalence > 0?'", 
     mar = c(2,0,3,2))
# this is higher than I think it should be .... but perhaps this is a function of our n ....
# or a function of the hypothesis test ... "Can we detect if prev > 0 | prev > 0?"
# which really only makes sense for k13

wald_ci <- function(n, p, alpha){
  z <- qnorm(alpha)
  #plot(sqrt(p * (1 - p) / n))
  pm <- z * sqrt(p * (1 - p) / n)
  return(list(upper = p + pm,
              lower = p - pm))
}

wilson_ci <- function(n, p, alpha){
  z <- qnorm(alpha)
  centre <- p + z**2/(2*n)
  pm <- (z / (2*n) * sqrt(4*n*p*(1 - p) + z**2))
  return(list(upper = (1 + z**2/n)**(-1)*(p + z**2/(2*n) + pm),
              lower = (1 + z**2/n)**(-1)*(p + z**2/(2*n) - pm)))
}

tmp <- wald_ci(n, p, 0.9)
# low sample size and low p mean that this struggles ...

tmp <- wilson_ci(n, p, 0.95)

par(mfrow = c(1,3))
plot(tmp$lower, main = "Lower")
plot(tmp$upper, main = "Upper")
# this is what I'm trying to get to!
plot(tmp$upper - tmp$lower, main = "Confidence interval width")
# looks quite a bit like 1/sqrt(n), so p isn't doing too much?
# would still like to get to the idea that we know *less* where our estimated p is low?
plot(1/sqrt(n))
# confidence interval width is mostly constant! duh! it's in the graph you just looked at!

# no idea what's going on down here I'm cooked but can't trust outputs until I figure it out
xxx <- seq(-8, 8, length.out = 100)
lll <- (1 + exp(-xxx))**(-1)

pal = c("darkblue", "mediumblue", "royalblue1")
{plot(xxx, lll, type = "l", 
     ylim = c(-0.1, 1.1), xlim = c(-7, 7),
     bty="n", xlab = "", ylab = "p", col=pal[1])
  abline(h = c(0,1))
  abline(h = 0.5, lty=2)
  lines(rep(0, 2), c(0, 1))
  # axis(2, line = -21.2) # no doubt this will move .. could just put at far left
  # mtext("p", 2, line = -19)
  tmp = wald_ci(50, lll, 0.95)
  lines(xxx, tmp$lower, col = pal[2])
  lines(xxx, tmp$upper, col = pal[2])
  text(1, 0.53, "n = 50", col = pal[2])
  
  tmp = wald_ci(5, lll, 0.95)
  lines(xxx, tmp$lower, col = pal[3])
  lines(xxx, tmp$upper, col = pal[3])
  text(1.5, 0.43, "n = 5", col = pal[3])
  mtext("Wald CI", adj = 0.25, font = 2, cex = 1.1)
  text(6, 0.85, "narrow CI\nat p ~ 0 or 1", col = "red3")
  text(-0.5, 1.1, "overshoot when n is small\nand p ~ 0 or 1", col = "red3")}


{plot(xxx, lll, type = "l", 
      ylim = c(-0.1, 1.1), xlim = c(-7, 7),
      bty="n", xlab = "", ylab = "p", col=pal[1])
  abline(h = c(0,1))
  abline(h = 0.5, lty=2)
  lines(rep(0, 2), c(0, 1))
  # axis(2, line = -21.2) # no doubt this will move .. could just put at far left
  # mtext("p", 2, line = -19)
  tmp = wilson_ci(50, lll, 0.95)
  lines(xxx, tmp$lower, col = pal[2])
  lines(xxx, tmp$upper, col = pal[2])
  text(1.1, 0.53, "n = 50", col = pal[2])
  
  tmp = wilson_ci(5, lll, 0.95)
  lines(xxx, tmp$lower, col = pal[3])
  lines(xxx, tmp$upper, col = pal[3])
  text(1.8, 0.43, "n = 5", col = pal[3])
  mtext("Wilson CI", adj = 0.25, font = 2, cex = 1.1)
  text(6, 0.7, "CI width ~ n", col = "red3")
  text(6, 1.1, "No overshoot", col = "red3")}














