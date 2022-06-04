
iter = function(Yobs, X, K, nchain,
                tau, nu0, sigma0, alpha1p, p1o, alpha0p, p0o, alpha1, alpha2,
                GAMMA_vec, p1_vec, p0_vec, Theta_vec, BETA_vec, sigma_vec, Y_vec,
                cur, update_BETA, dependency, strength) {
  for (ich in 1:nchain) {
    Y = Y_vec[[ich]][[cur]]
    p1 = p1_vec[[ich]][[cur]]
    p0 = p0_vec[[ich]][[cur]]
    GAMMA = GAMMA_vec[[ich]][[cur]]
    Theta = Theta_vec[[ich]][[cur]]
    sigma = sigma_vec[[ich]][cur]
    BETA = BETA_vec[[ich]][[length(BETA_vec[[ich]])]]

    l1 = lambda1(p1)
    l2 = lambda2(p0)
    Gammavec = GammaGprod(Zprod(GAMMA, l1, l2, Theta),
                          GAMMA, l1, l2, Theta)
    for (k in 1:K) {
      locur = 1
      repeat{
        GAMMAnew = Gamma_prop(GAMMA, k)
        kres = pGamma_ratio(GAMMAnew, GAMMA, k, X[[k]], Y[[k]],
                            tau[k], sigma0, nu0, Gammavec,
                            l1, l2, Theta)
        locur = locur + 1
        alpha_ratio = runif(1, 0.0001, 0.9999)
        if(kres$ratio > alpha_ratio) {
          GAMMA = GAMMAnew
          Gammavec = kres$Gammavec
          break
        } else if(locur > 10) {
          break
        }
      }
    }
    GAMMA_vec[[ich]][[cur + 1]] = GAMMA

    if (dependency) {
      locur = 1
      repeat{
        p1new = P1_prop(p1, alpha1p)
        l1 = lambda1(p1new)
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2, Theta), GAMMA, l1, l2, Theta)
        locur = locur + 1
        new_ratio = pP1_ratio(p1new, p1, alpha1p, p1o[[ich]], Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0.0001, 0.9999)
        if(new_ratio > alpha_ratio) {
          p1 = p1new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      p1_vec[[ich]][[cur + 1]] = p1

      locur = 1
      repeat{
        p0new = P0_prop(p0, alpha0p)
        l2 = lambda2(p0new)
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2, Theta), GAMMA, l1, l2, Theta)
        locur = locur + 1
        new_ratio = pP0_ratio(p0new, p0, alpha0p, p0o[[ich]], Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0.0001, 0.9999)
        if(new_ratio > alpha_ratio) {
          p0 = p0new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      p0_vec[[ich]][[cur + 1]] = p0
    } else {
      locur = 1
      repeat{
        p0new = P0_prop(p0, alpha0p)
        p1new = p0new
        l1 = lambda1(p1new)
        l2 = lambda2(p0new)
        Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2, Theta), GAMMA, l1, l2, Theta)
        locur = locur + 1
        new_ratio = pP0_ratio(p0new, p0, alpha0p, p0o[[ich]], Gammavec_new, Gammavec)
        alpha_ratio = runif(1, 0.0001, 0.9999)
        if(new_ratio > alpha_ratio) {
          p0 = p0new
          Gammavec = Gammavec_new
          break
        } else if(locur > 10) {
          break
        }
      }
      p0_vec[[ich]][[cur + 1]] = p0
      p1_vec[[ich]][[cur + 1]] = p1
    }

    if (strength) {
      for (k in 1:(K - 1)) {
        for (kp in (k + 1):K) {
          Thetanew = Theta_prop(Theta, alpha1, k, kp)
          if (dependency) {
            Gammavec_new = GammaGprod(Zprod(GAMMA, l1, l2, Thetanew), GAMMA, l1, l2, Thetanew)
          } else {
            Gammavec_new = GammaGprod2(Zprod2(GAMMA, l1, Thetanew), GAMMA, l1, Thetanew)
          }
          locur = locur + 1
          new_ratio = pTheta_ratio(Thetanew, Theta, GAMMA, k, kp, alpha1, alpha2, Gammavec_new, Gammavec)
          alpha_ratio = runif(1, 0.0001, 0.9999)
          if(new_ratio > alpha_ratio) {
            Theta = Thetanew
            Gammavec = Gammavec_new
          }
        }
      }
    } else {
      Theta = matrix(0, K, K)
    }
    Theta_vec[[ich]][[cur + 1]] = Theta

    sigma = sigma_prop(sigma, GAMMA, X, Y, tau, nu0, sigma0)
    sigma_vec[[ich]] = c(sigma_vec[[ich]], sigma)

    for (k in 1:K) {
      Y = Yk_prop(Y, Yobs, X, GAMMA, tau, nu0, sigma0, k)
    }
    Y_vec[[ich]][[cur + 1]] = Y
  }

  cur = cur + 1
  result = list(cur = cur,
                GAMMA_vec = GAMMA_vec,
                p1_vec = p1_vec,
                p0_vec = p0_vec,
                Theta_vec = Theta_vec,
                BETA_vec = BETA_vec,
                sigma_vec = sigma_vec,
                Y_vec = Y_vec)
  return(result)
}


#' Main model fitting function
#'
#' This is the main model fitting function for Bayesian framework with ISP and MRF priors.
#' @details
#' This function is used to fit a Baeysian regression model with ISP and MRF priors.
#'
#' @param Y A list of response values in regressions
#' @param X A list of design matrices in regressions
#' @param iteration The number of MCMC iterations
#' @param burnin The number of burn-in iterations
#' @param nchain The number of MCMC chains
#' @param dependency Whether to use the ISP prior; A value of TRUE means creating dependencies among SNPs
#' @param strength Whether to borrow strength across regressions using MRF prior
#' @return A large list of all the parameter values
#' @examples
#' library(mvtnorm)
#' K = 3
#' P = 500
#' Nk = c(50, 100, 150)
#' GAMMA = matrix(0, nrow = K, ncol = P)
#' GAMMA[1, 1:8] = 1
#' GAMMA[2, 1:8] = 1
#' GAMMA[3, 1:8] = 1
#' Y = lapply(Nk, function(l) rep(0, l))
#' X = lapply(Nk, function(l) matrix(0, l, P))
#' sigma = lapply(Nk, function(l) rep(0, l))
#' SigmaX = 0.05 ^ abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
#' BETA = matrix(0, nrow = K, ncol = P)
#' sigma = lapply(Nk, function(l) rep(0, l))
#'
#' for (k in 1:K) {
#'  X[[k]] = matrix(rmvnorm(Nk[k], sigma = SigmaX), nrow = Nk[k], byrow = TRUE)
#'  BETA[k, which(GAMMA[k, ] == 1)] = c(2, -2, 2, -2, 2, -2, 2, -2)
#'  sigma[[k]] = rnorm(Nk[k], 0, 1)
#'  Y[[k]] = X[[k]] %*% (GAMMA[k, ] * BETA[k, ]) + sigma[[k]]
#'  Y[[k]] = as.numeric(Y[[k]] > 0)
#' }
#'
#' t0 = Sys.time()
#' res = Main(Y, X, iteration = 5000, burnin = 3500, nchain = 2)
#' Sys.time() - t0
#'
#' @export

Main = function(Y, X, iteration, burnin, nchain,
                tau = NULL, nu0 = 5, sigma0 = NULL,
                alpha1p = 2, p1o = NULL, alpha0p = 2, p0o = NULL,
                alpha1 = NULL, alpha2 = NULL,
                dependency = TRUE, strength = TRUE) {
  K = length(Y)
  G = dim(X[[1]])[2]
  if(is.null(tau)) {
    tau = unlist(lapply(Y, function(y) sd(y)))
  }
  if(is.null(sigma0)) {
    sigma0 = sd(unlist(Y)) ^ 2
  }
  if(is.null(p1o)) {
    p1o = rep(0.25, G)
  }
  if(is.null(p0o)) {
    p0o = rep(0.1, G)
  }
  if(is.null(alpha1)) {
    alpha1 = 5
  }
  if(is.null(alpha2)) {
    alpha2 = 1
  }

  p0_vec = list()
  for (ich in 1:nchain) {
    p0_vec[[ich]] = list()
    if (length(unique(p0o)) == 1) {
      p0_vec[[ich]][[1]] = rep(rbeta(1, alpha0p, alpha0p * (1 - p0o[1]) / p0o[1]), G)
    } else {
      p0_vec[[ich]][[1]] = unlist(lapply(p0o, function(x) rbeta(1, alpha0p, alpha0p * (1 - x) / x)))
    }
  }
  if (dependency) {
    p1_vec = list()
    for (ich in 1:nchain) {
      p1_vec[[ich]] = list()
      if (length(unique(p1o)) == 1) {
        p1_vec[[ich]][[1]] = rep(rbeta(1, alpha1p, alpha1p * (1 - p1o[1]) / p1o[1]), G)
      } else {
        p1_vec[[ich]][[1]] = unlist(lapply(p1o, function(x) rbeta(1, alpha1p, alpha1p * (1 - x) / x)))
      }
    }
  } else {
    p1_vec = p0_vec
  }
  Theta_vec = list()
  for (ich in 1:nchain) {
    Theta_vec[[ich]] = list()
    if (strength) {
      Theta_vec[[ich]][[1]] = matrix(runif(1, 0, 1.5), K, K)
      diag(Theta_vec[[ich]][[1]]) = 0
    } else {
      Theta_vec[[ich]][[1]] = matrix(0, K, K)
    }
  }
  GAMMA_vec = list()
  for (ich in 1:nchain) {
    GAMMA_vec[[ich]] = list()
    GAMMA = matrix(0, K, G)
    GAMMA[, 1] = sample(c(0, 1), K, replace = TRUE, prob = c(1 - p1_vec[[ich]][[1]][1], p1_vec[[ich]][[1]][1]))
    for (g in 2:G) {
      for (k in 1:K) {
        if(GAMMA[k, g - 1] == 1) {
          GAMMA[k, g] = sample(c(0, 1), 1, prob = c(1 - p1_vec[[ich]][[1]][g], p1_vec[[ich]][[1]][g]))
        } else {
          GAMMA[k, g] = sample(c(0, 1), 1, prob = c(1 - p0_vec[[ich]][[1]][g], p0_vec[[ich]][[1]][g]))
        }
      }
    }
    GAMMA_vec[[ich]][[1]] = GAMMA
  }
  BETA_vec = list()
  sigma_vec = list()
  for (ich in 1:nchain) {
    BETA_vec[[ich]] = list()
    BETA_vec[[ich]][[1]] = matrix(0, K, G)
    sigma_vec[[ich]] = vector()
    sigma_vec[[ich]][1] = sigma_prop(sigma0, GAMMA_vec[[ich]][[1]], X, Y, tau, nu0, sigma0)
  }
  l1 = list()
  l2 = list()
  for (ich in 1:nchain) {
    l1[[ich]] = as.numeric(lambda1(p1_vec[[ich]][[1]]))
    l2[[ich]] = as.numeric(lambda2(p0_vec[[ich]][[1]]))
  }
  Gammavec = list()
  for (ich in 1:nchain) {
    if (dependency) {
      Gammavec[[ich]] = GammaGprod(Zprod(GAMMA_vec[[ich]][[1]], l1[[ich]], l2[[ich]], Theta_vec[[ich]][[1]]),
                                   GAMMA_vec[[ich]][[1]], l1[[ich]], l2[[ich]], Theta_vec[[ich]][[1]])
    } else {
      Gammavec[[ich]] = GammaGprod2(Zprod2(GAMMA_vec[[ich]][[1]], l1[[ich]], Theta_vec[[ich]][[1]]),
                                    GAMMA_vec[[ich]][[1]], l1[[ich]], Theta_vec[[ich]][[1]])
    }
  }
  Y_vec = list()
  for (ich in 1:nchain) {
    Y_vec[[ich]] = list()
    Y_vec[[ich]][[1]] = Y
  }

  cur = 1
  for (index in 1:iteration) {
    if(cur <= burnin) {
      update_BETA = FALSE
    } else if(cur > burnin)  {
      update_BETA = TRUE
    }
    iter_res = iter(Yobs = Y, X, K, nchain,
                    tau, nu0, sigma0, alpha1p, p1o, alpha0p, p0o, alpha1, alpha2,
                    GAMMA_vec, p1_vec, p0_vec, Theta_vec, BETA_vec, sigma_vec, Y_vec,
                    cur, update_BETA, dependency, strength)
    cur = iter_res$cur
    GAMMA_vec = iter_res$GAMMA_vec
    p1_vec = iter_res$p1_vec
    p0_vec = iter_res$p0_vec
    Theta_vec = iter_res$Theta_vec
    BETA_vec = iter_res$BETA_vec
    sigma_vec = iter_res$sigma_vec
    Y_vec = iter_res$Y_vec
  }

  return(iter_res)
}
