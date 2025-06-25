pwr_binom <- function(n, p){
  # probability k > 0 given n, p
  1 - (1 - p)^n
}

p <- preds$`2019_post_mean`

n <- test_dens_masked %>%
  project(p) %>% # for all years ...
  mask(p)

p[cells(p)] <- 0.1
n[n < 1] <- NA

pwr <- pwr_binom(n, p)

plot(sqrt(n), main="Surv intensity")
plot(pwr, main="Power")
# this is higher than I think it should be .... but perhaps this is a function of our n ....