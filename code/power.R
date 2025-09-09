preds <- rast("output/k13_marcse/gneiting_ahmc/preds_all.tif")
ras <- rast("output/k13_marcse/surveillance_effort_k13_marcse.grd")

pwr_binom <- function(n, p){
  # probability k > 0 given n, p
  1 - (1 - p)^n
}

p <- preds$`2022_post_median`

# from surveillance_effort
# gf <- focalMat(ras, 0.8, "Gauss")
# n <- aggregate(ras, fact = 4, sum, na.rm = TRUE) %>%
#   focal(gf, na.rm=TRUE) %>%
#   project(p) %>%
#   mask(p)

# aggregate makes more sense here ....
# n <- test_dens %>%
#   aggregate(fact = 4, sum, na.rm=TRUE) %>%
#   project(p) %>%
#   mask(p)

# aggregate, then disaggregate to adjust units to tests per sqkmish
n <- ras %>%
  aggregate(fact = 2, sum, na.rm = TRUE) %>%
  project(p) %>%
  mask(p)

plot(n)

# could also use extant test_dens + 1?

# shouldn't be a problem ...?
# n[n < 1] <- NA

pwr <- pwr_binom(n, p)

mar = c(4, 2, 3, 1)
{png("~/Desktop/presentations/marcse/power.png", 
     height = 800, width = 2400, pointsize = 40)
par(mfrow = c(1,3))
plot(trim(p), main="Estimated prevalence", mar = mar, plg = list(x = "bottom"))
plot(trim(n), main="Surv intensity", mar = mar, plg = list(x = "bottom"))
plot(trim(pwr), main='Statistical power:\n "Pr(detect prevalence > 0 | prevalence > 0)"', 
     mar = mar, plg = list(x = "bottom"))
dev.off()}
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

wilson_ci <- function(n, p, alpha, epsilon = 1e-5){
  n[n == 0] = epsilon
  z <- qnorm(alpha)
  centre <- p + (z**2)/(2*n)
  scale <- (1 + (z**2)/n)**(-1)
  
  pm <- (z / (2*n)) * sqrt(4*n*p*(1 - p) + z**2)
  
  # oops had some brackets in the wrong place here >:)
  return(list(upper = scale*(centre + pm),
              lower = scale*(centre - pm)))
  
  # end up with -ves being given to sqrt?
}

n <- n + 0.1
tmp <- wald_ci(n, p, 0.95)
# low sample size and low p mean that this struggles ...

tmp <- wilson_ci(n, p, 0.95)

n[n == 0] <- 1e-5
df <- xyFromCell(p, cell = cells(p))
ras <- rast(list(rast(tmp), p, n))
ras$`CI width` <- ras$upper - ras$lower
# ras$lower[ras$lower <= 0] = NA
n = n + 1
ras$inv_eff <- 1/sqrt(n)

df <- cbind(df, terra::extract(ras, df))
names(df) <- c("Longitude", "Latitude", "CI upper", "CI lower", 
               "Estimated prevalence", "surveillance intensity", "CI width", "inv_eff")
df <- df %>%
  pivot_longer(cols = c("CI upper", "CI lower", 
                        "Estimated prevalence", "surveillance intensity", 
                        "CI width", "inv_eff"),
               names_to = "lyr")

p1 <- ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = Longitude, y = Latitude, fill = value), 
            data = df %>% filter(!lyr %in% c("surveillance intensity", "inv_eff"))) +
  scale_fill_viridis_c(name = "Prevalence") +
  facet_wrap(~lyr)
p2 <- ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = Longitude, y = Latitude, fill = value), 
            data = df %>% filter(lyr == "surveillance intensity")) +
  scale_fill_viridis_c(name = "Tests per\n100sqkm") +
  labs(title = "Surveillance intensity")
p3 <- ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = Longitude, y = Latitude, fill = value), 
            data = df %>% filter(lyr == "inv_eff")) +
  scale_fill_viridis_c(name = "Inverse effort") +
  labs(title = "1/sqrt(surveillance intensity)")
p4 <- plot_grid(p2, p3, ncol = 1)
plot_grid(p1, p4, rel_widths = c(1, 0.62))
ggsave("~/Desktop/presentations/marcse/ci_wald.png", height = 4, width = 6, scale = 2)


{par(mfrow = c(2,2))
plot(tmp$lower, main = "Lower")
plot(tmp$upper, main = "Upper")
# this is what I'm trying to get to!
plot(tmp$upper - tmp$lower, main = "Confidence interval width")
# looks quite a bit like 1/sqrt(n), so p isn't doing too much?
# would still like to get to the idea that we know *less* where our estimated p is low?
plot(1/sqrt(n), main = "1/sqrt(n)")}
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


plot(n)
tmp <- n
tmp[tmp < 0.05] <- NA
plot(tmp)


# overdispersion:
# V(Y_i) = \sigma^2 * \mu_i(n_i - \mu_i) / n_i
# \hat\sigma**2 = X**2 / (N - p), 
# X**2 is Pearson goodness of fit, 
# p == number of params, N == number of records
# e.g.....
8600 / (680 - p)
# I would say that's overdispersed ??
# If I've done this right?








