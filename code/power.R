pwr_binom <- function(n, p){
  # probability k > 0 given n, p
  1 - (1 - p)^n
}

p <- preds$`2022_post_mean`

# from surveillance_effort
gf <- focalMat(ras, 1, "Gauss")
n <- aggregate(ras, fact = 4, sum, na.rm = TRUE) %>%
  focal(gf, na.rm=TRUE) %>%
  project(p) %>%
  mask(p)

# aggregate makes more sense here ....
n <- test_dens %>%
  aggregate(fact = 4, sum, na.rm=TRUE) %>%
  project(p) %>%
  mask(p)
n <- n + 1

# could also use extant test_dens + 1?

n[n < 1] <- NA

pwr <- pwr_binom(n, p)

plot(sqrt(n), main="Surv intensity")
plot(pwr, main="Power")
# this is higher than I think it should be .... but perhaps this is a function of our n ....

wald_ci <- function(n, p, alpha){
  z <- qnorm(alpha)
  plot(sqrt(p * (1 - p) / n))
  p + c(1, -1) * z * sqrt(p * (1 - p) / n)
}

wilson_ci <- function(n, p, alpha){
  z <- qnorm(alpha)
  return((1 + z**2/n)**(-1)*(p + z**2/(2*n) + c(-1, 1)*(z / (2*n) * sqrt(4*n*p*(1 - p) + z**2))))
}

tmp <- wald_ci(n, p, 0.9)
# low sample size and low p mean that this struggles ...

tmp <- wilson_ci(n, p, 0.95)

par(mfrow = c(1,3))
plot(tmp$lyr1, main = "Lower")
plot(tmp$lyr2, main = "Upper")
# this is what I'm trying to get to!
plot(tmp$lyr2 - tmp$lyr1, main = "Confidence interval width")
# looks quite a bit like 1/sqrt(n), so p isn't doing too much?
# would still like to get to the idea that we know *less* where our estimated p is low?
plot(1/sqrt(n))
# confidence interval width is mostly constant! duh! it's in the graph you just looked at!

# no idea what's going on down here I'm cooked but can't trust outputs until I figure it out
ppp <- seq(0,1, length.out = 100)
tmp <- wilson_ci(100, ppp, 0.95)
plot(ppp, tmp[seq])







