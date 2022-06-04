
lambda1 = function(p1) {
  unlist(lapply(p1, function(x) log(x) - log(1 - x)))
}


lambda2 = function(p0) {
  unlist(lapply(p0, function(x) log(x) - log(1 - x)))
}


expGamma = function(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta) {
  sum = lambda1g * gamma_g %*% gamma_g_1
  sum = sum + lambda2g * gamma_g %*% (1 - gamma_g_1)
  sum = sum + gamma_g %*% Theta %*% gamma_g
  return(exp(sum))
}


expGamma1 = function(gamma_g, lambda1g, Theta) {
  sum = lambda1g * gamma_g %*% gamma_g
  sum = sum + gamma_g %*% Theta %*% gamma_g
  return(exp(sum))
}


Zloop = function(vecs, gamma_g_1, lambda1g, lambda2g, Theta) {
  sum = apply(vecs, 2, function(x) expGamma(x, gamma_g_1, lambda1g, lambda2g, Theta))
  return(sum(sum))
}


Zloop1 = function(vecs, lambda1g, Theta) {
  sum = apply(vecs, 2, function(x) expGamma1(as.numeric(x), lambda1g, Theta))
  return(sum(sum))
}


Z = function(K, gamma_g_1, lambda1g, lambda2g, Theta) {
  vecs = t(as.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE))))
  sum = Zloop(vecs, gamma_g_1, lambda1g, lambda2g, Theta)
  return(sum)
}


pGamma_g = function(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta) {
  K = length(gamma_g)
  ratio = expGamma(gamma_g, gamma_g_1, lambda1g, lambda2g, Theta) / Z(K, gamma_g_1, lambda1g, lambda2g, Theta)
  return(ratio)
}


Z1 = function(K, lambda1g, Theta) {
  vecs = t(as.matrix(expand.grid(replicate(K, 0:1, simplify = FALSE))))
  sum = Zloop1(vecs, lambda1g, Theta)
  return(sum)
}


pGamma_1 = function(gamma_g, lambda1g, Theta) {
  K = length(gamma_g)
  ratio = expGamma1(gamma_g, lambda1g, Theta) / Z1(K, lambda1g, Theta)
  return(ratio)
}


Zprod = function(GAMMA, l1, l2, Theta) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Zvec = c(Z1(K, l1[1], Theta),
           unlist(lapply(2:G, function(g) Z(K, GAMMA[, g - 1], l1[g], l2[g], Theta) )))
  return(Zvec)
}


Zprod2 = function(GAMMA, l1, Theta) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Zvec = unlist(lapply(1:G, function(g) Z(K, GAMMA[, g], l1[g], 0, Theta) ))
  return(Zvec)
}


GammaGprod = function(Zvec, GAMMA, l1, l2, Theta) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Gammavec = c(expGamma1(GAMMA[, 1], l1[1], Theta),
               unlist(lapply(2:G, function(g) expGamma(GAMMA[, g], GAMMA[, g - 1], l1[g], l2[g], Theta) )))
  Gammavec = exp(log(Gammavec) - log(Zvec))
  return(Gammavec)
}


GammaGprod2 = function(Zvec, GAMMA, l1, Theta) {
  G = dim(GAMMA)[2]
  K = dim(GAMMA)[1]

  Gammavec = unlist(lapply(1:G, function(g) expGamma(GAMMA[, g], GAMMA[, g], l1[g], 0, Theta) ))
  Gammavec = exp(log(Gammavec) - log(Zvec))
  return(Gammavec)
}


pGamma_ratio = function(GAMMAnew, GAMMAold, k, Xk, Yk, tauk, sigma0, nu0, Gammavec_old,
                        l1, l2, Theta) {
  Gammavec_new = Gammavec_old
  c = GAMMAnew[k, ] - GAMMAold[k, ]
  c = which(c != 0)
  for (i in 1:length(c)) {
    if (c[i] == 1) {
      Gammavec_new[1] = pGamma_1(GAMMAnew[, 1], l1[1], Theta)
    } else if (c[i] < length(l1)) {
      Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i] - 1],
                                    l1[c[i]], l2[c[i]], Theta)
      Gammavec_new[c[i] + 1] = pGamma_g(GAMMAnew[, c[i] + 1], GAMMAnew[, c[i]],
                                        l1[c[i] + 1], l2[c[i] + 1], Theta)
    } else {
      Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i] - 1],
                                    l1[c[i]], l2[c[i]], Theta)
    }
  }

  gamma_new = which(GAMMAnew[k, ] == 1)
  gamma_old = which(GAMMAold[k, ] == 1)
  Sigma_gammanew = tauk * Xk[, gamma_new] %*% t(Xk[, gamma_new])
  Sigma_gammaold = tauk * Xk[, gamma_old] %*% t(Xk[, gamma_old])
  diag(Sigma_gammanew) = diag(Sigma_gammanew) + 1
  diag(Sigma_gammaold) = diag(Sigma_gammaold) + 1
  Sigma_gammanew = sigma0 * Sigma_gammanew
  Sigma_gammaold = sigma0 * Sigma_gammaold
  mulTnew = dmvt(as.numeric(Yk), sigma = Sigma_gammanew, df = nu0, log = FALSE)
  mulTold = dmvt(as.numeric(Yk), sigma = Sigma_gammaold, df = nu0, log = FALSE)

  ratio = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * mulTnew / mulTold
  if(is.na(ratio)) ratio = 0
  return(list(ratio = ratio,
              Gammavec = Gammavec_new))
}


pGamma_ratio2 = function(GAMMAnew, GAMMAold, k, Xk, Yk, tauk, sigma0, nu0, Gammavec_old,
                         l1, Theta) {
  Gammavec_new = Gammavec_old
  c = GAMMAnew[k, ] - GAMMAold[k, ]
  c = which(c != 0)
  for (i in 1:length(c)) {
    Gammavec_new[c[i]] = pGamma_g(GAMMAnew[, c[i]], GAMMAnew[, c[i]],
                                  l1[c[i]], 0, Theta)
  }

  gamma_new = which(GAMMAnew[k, ] == 1)
  gamma_old = which(GAMMAold[k, ] == 1)
  Sigma_gammanew = tauk * Xk[, gamma_new] %*% t(Xk[, gamma_new])
  Sigma_gammaold = tauk * Xk[, gamma_old] %*% t(Xk[, gamma_old])
  diag(Sigma_gammanew) = diag(Sigma_gammanew) + 1
  diag(Sigma_gammaold) = diag(Sigma_gammaold) + 1
  Sigma_gammanew = sigma0 * Sigma_gammanew
  Sigma_gammaold = sigma0 * Sigma_gammaold
  mulTnew = dmvt(as.numeric(Yk), sigma = Sigma_gammanew, df = 2 * nu0, log = FALSE)
  mulTold = dmvt(as.numeric(Yk), sigma = Sigma_gammaold, df = 2 * nu0, log = FALSE)

  ratio = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * mulTnew / mulTold
  if(is.na(ratio)) ratio = 0
  return(list(ratio = ratio,
              Gammavec = Gammavec_new))
}


pTheta_ratio = function(Thetanew, Thetaold, GAMMA, k, kp, alpha1, alpha2, Gammavec_new, Gammavec_old) {
  thetanew = Thetanew[k, kp]
  thetaold = Thetaold[k, kp]
  res = dgamma(thetanew, shape = alpha1, rate = alpha2) /
    dgamma(thetaold, shape = alpha1, rate = alpha2)
  res = res * exp(sum(log(Gammavec_new) - log(Gammavec_old)))
  if(is.na(res)) res = 0
  return(res)
}


pP1_ratio = function(p1_new, p1_old, alpha1p, p1o, Gammavec_new, Gammavec_old) {
  if (length(unique(p1_new)) == 1) {
    dnew = dbeta(p1_new[1], alpha1p, alpha1p * (1 - p1o[1]) / p1o[1])
    dold = dbeta(p1_old[1], alpha1p, alpha1p * (1 - p1o[1]) / p1o[1])
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * dnew / dold
    if(is.na(res)) res = 0
    return(res)
  } else {
    dnew = rep(0, length(p1_new))
    dold = rep(0, length(p1_old))
    for (g in 1:length(p1_new)) {
      dnew[g] = dbeta(p1_new[g], alpha1p, alpha1p * (1 - p1o[g]) / p1o[g])
      dold[g] = dbeta(p1_old[g], alpha1p, alpha1p * (1 - p1o[g]) / p1o[g])
    }
    res = exp(sum(log(Gammavec_new * dnew) - log(Gammavec_old * dold)))
    if(is.na(res)) res = 0
    return(res)
  }
}


pP0_ratio = function(p0_new, p0_old, alpha0p, p0o, Gammavec_new, Gammavec_old) {
  if (length(unique(p0_new)) == 1) {
    dnew = dbeta(p0_new[1], alpha0p, alpha0p * (1 - p0o[1]) / p0o[1])
    dold = dbeta(p0_old[1], alpha0p, alpha0p * (1 - p0o[1]) / p0o[1])
    res = exp(sum(log(Gammavec_new) - log(Gammavec_old))) * dnew / dold
    if(is.na(res)) res = 0
    return(res)
  } else {
    dnew = rep(0, length(p0_new))
    dold = rep(0, length(p0_old))
    for (g in 1:length(p0_new)) {
      dnew[g] = dbeta(p0_new[g], alpha0p, alpha0p * (1 - p0o[g]) / p0o[g])
      dold[g] = dbeta(p0_old[g], alpha0p, alpha0p * (1 - p0o[g]) / p0o[g])
    }
    res = exp(sum(log(Gammavec_new * dnew) - log(Gammavec_old * dold)))
    if(is.na(res)) res = 0
    return(res)
  }
}


Gamma_prop = function(GAMMA, k) {
  Gamma = GAMMA[k, ]
  sors = sample(c(0, 1), 1)
  if (sors == 0 & length(which(Gamma == 1)) > 8) {
    ones = which(GAMMA[k, ] == 1)
    zeros = which(GAMMA[k, ] == 0)
    Gamma[sample(ones, 1)] = 0
    Gamma[sample(zeros, 1)] = 1
  } else if (sors == 1 & length(which(Gamma == 1)) > 8) {
    c1 = sample(1:length(Gamma), 1)
    Gamma[c1] = 1 - Gamma[c1]
  } else if (length(which(Gamma == 1)) <= 8) {
    zeros = which(GAMMA[k, ] == 0)
    Gamma[sample(zeros, 1)] = 1
  }

  GAMMA[k, ] = Gamma
  return(GAMMA)
}


Theta_prop = function(Theta, alpha1, k, kp) {
  theta = Theta[k, kp]
  theta = rgamma(1, shape = alpha1, rate = alpha1 / theta)
  Theta[k, kp] = theta
  Theta[kp, k] = theta
  return(Theta)
}


P1_prop = function(p1, alpha1p) {
  if (length(unique(p1)) == 1) {
    p1new = rbeta(1, alpha1p, alpha1p * (1 - p1[1]) / p1[1])
    return(rep(p1new, length(p1)))
  } else {
    p1new = rep(0, length(p1))
    for (i in 1:length(p1)) p1new[i] = rbeta(1, alpha1p, alpha1p * (1 - p1[i]) / p1[i])
    return(p1new)
  }
}


P0_prop = function(p0, alpha0p) {
  if (length(unique(p0)) == 1) {
    p0new = rbeta(1, alpha0p, alpha0p * (1 - p0[1]) / p0[1])
    return(rep(p0new, length(p0)))
  } else {
    p0new = rep(0, length(p0))
    for (i in 1:length(p0)) p0new[i] = rbeta(1, alpha0p, alpha0p * (1 - p0[i]) / p0[i])
    return(p0new)
  }
}


sigma_prop = function(sigma, GAMMA, X, Y, tau, nu0, sigma0) {
  shape = (sum(unlist(lapply(Y, length))) + nu0) / 2
  scale = 0
  for (k in 1:dim(GAMMA)[1]) {
    Xk = X[[k]]
    Yk = Y[[k]]
    Gamma = which(GAMMA[k, ] == 1)
    SigmaY = tau[k] * Xk[, Gamma] %*% t(Xk[, Gamma]) + diag(length(Yk))
    scale = scale + Yk %*% solve(SigmaY) %*% Yk
  }
  scale = scale / 2 + nu0 * sigma0 / 2
  return(MCMCpack::rinvgamma(1, shape = shape, scale = scale))
}


Yk_prop = function(Y, Yobs, X, GAMMA, tau, nu0, sigma0, k) {
  Xk = X[[k]]
  Yk = Y[[k]]
  yobs = Yobs[[k]]
  gamma = which(GAMMA[k, ] == 1)
  Sigma_gamma = tau[k] * Xk[, gamma] %*% t(Xk[, gamma])
  diag(Sigma_gamma) = diag(Sigma_gamma) + 1
  Sigma_gamma = sigma0 * Sigma_gamma
  ll = rep(0, nrow(Sigma_gamma))
  ll[yobs == 0] = -Inf
  uu = rep(0, nrow(Sigma_gamma))
  uu[yobs == 1] = Inf
  Yk = rtmvt(1, mean = rep(0, nrow(Sigma_gamma)), sigma = Sigma_gamma,
             df = nu0, lower = ll, upper = uu, algorithm = "gibbs")
  Y[[k]] = as.numeric(Yk)
  return(Y)
}

