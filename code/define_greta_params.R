#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
define_greta_parameters <- function(n_snp,
                                    n_latent) {
  
  
  beta_mu <- define_greta_beta_mu()
  beta_sd <- define_greta_beta_sd()
  
  list(
    # individual probabilities of failure per SNP
    q = define_greta_q(n_snp = n_snp),
    # SNP - latent factor loadings
    loadings = define_greta_lambda(n_snp = n_snp,
                                   n_latent = n_latent),
    # hierarchical regression coefficients
    beta_mu = beta_mu,
    beta_sd = beta_sd,
    beta = define_greta_beta(n_snp = n_snp,
                             mu = beta_mu,
                             sd = beta_sd)
  )
  
}


#' Create a greta array for the hierarchical standard deviation of the
#' regression coefficients between SNPs and covariates
#'
#' @title Define 'sd' Parameter
#' @return a greta array of dimension 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta_sd <- function() {
  
  greta::normal(0, 1, dim = 3, truncation = c(0, Inf))
  
}


#' Create a greta array for the hierarchical mean of the regression coefficients
#' between SNPs and covariates
#'
#' @title Define 'mu' Parameter
#' @return a greta array of dimension 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta_mu <- function() {
  
  greta::normal(0, 1, dim = 3)
  
}


#' Create a greta array for the hierarchical regression coefficients between
#' SNPs and covariates
#'
#' @title Define 'beta' Parameter
#' @param n_snp number of SNPs
#' @return a greta array of dimension `n_snp` x 3
#' @importFrom greta normal
#' @author Nick Golding
define_greta_beta <- function(n_snp,
                              mu,
                              sd) {
  
  # use hierarchical decentring specification, to remove prior correlation from
  # the posterior
  beta_raw <- greta::normal(0, 1, dim = c(n_snp, 3))
  beta_raw_sd <- greta::sweep(beta_raw, 2, sd, FUN = "*")
  beta <- greta::sweep(beta_raw_sd, 2, mu, FUN = "+")
  beta
  
}


#' Create a greta array for the loadings between SNPs and latent factors
#'
#' @title Define 'Lambda' Parameter
#' @param n_snp number of SNPs
#' @param n_latent number of latent factors
#' @return a greta array of dimension `n_snp` x `n_latent`
#' @importFrom greta zeros normal
#' @author Nick Golding
define_greta_lambda <- function(n_snp = n_snp, n_latent = n_latent) {
  
  # create lambda, subject to constraints to prevent rotational invariance when
  # fitting
  # upper triangular zero
  loadings <- greta::zeros(n_snp, n_latent)
  # diagonal positive
  diag(loadings) <- greta::normal(0, 1, truncation = c(0, Inf), dim = n_latent)
  # lower triangular unconstrained
  lower <- lower.tri(loadings)
  loadings[lower] <- greta::normal(0, 1, dim = sum(lower))
  
  loadings
}

#' Create a greta array for the individual probabilities of treatment failure
#' per SNP
#'
#' @title Define 'q' Parameter
#' @param n_snp number of SNPs
#' @return a greta array of dimension `n_snp`
#' @author Nick Golding
define_greta_q <- function(n_snp) {
  
  greta::cauchy(0, 0.1,
                dim = n_snp,
                truncation = c(0, 1))
  
}
