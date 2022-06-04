
#' Cross-validated prediction function
#'
#' This is the prediction function with built-in cross validation.
#' @details
#' This function performs cross-validation and generates predicted values.
#'
#' @param res The model output from the Main() function
#' @param Y A list of response values in regressions
#' @param X A list of design matrices in regressions
#' @param burnin The number of burn-in iterations
#' @param nfold The number of folds for cross validation
#' @param ich The index of chain used for prediction (e.g., ich = 1)
#' @return A list of predicted probabilities for each subject to be a case
#' @examples
#' Predict(res, Y, X, 3500)
#'
#' @export

Predict = function(res, Y, X, burnin, nfold = 5, ich = 1) {
  require(caret)
  require(mvtnorm)
  GAMMAres = res$GAMMA_vec[[ich]][(burnin + 1):length(res$GAMMA_vec[[ich]])]
  Yaug = list()
  for(k in 1:K) {
    Yres = lapply(res$Y_vec[[ich]], function(x) x[[k]])
    Yres = Yres[(burnin + 1):length(res$Y_vec[[ich]])]
    Yaug[[k]] = Reduce("+", Yres) / length(Yres)
  }
  Ypred = lapply(1:K, function(k) rep(0, length(Yaug[[k]])))
  sigmahat = mean(res$sigma_vec[[ich]][(burnin + 1):length(res$GAMMA_vec[[ich]])])
  folds = lapply(Y, function(x) createFolds(x, k = nfold, list = FALSE, returnTrain = FALSE))

  for (k in 1:K) {
    for (m in 1:nfold) {
      Xtrain = X[[k]][folds[[k]] != m, ]
      Ytrain = Y[[k]][folds[[k]] != m]
      Xtest = X[[k]][folds[[k]] == m, ]
      Yhtest = Yaug[[k]][folds[[k]] == m]
      W = rep(0, length(GAMMAres))

      pred_list = list()
      for (i in 1:length(GAMMAres)) {
        gamma = which(GAMMAres[[i]][k, ] == 1)
        fit = glmnet(Xtrain[, gamma], Ytrain, alpha = 1, family = "binomial", intercept = FALSE)
        Bhat = as.numeric(coef(fit, s = fit$lambda[length(fit$lambda)]))[-1]
        if (length(gamma) == 1) {
          pred_list[[i]] = Xtest[, gamma] * Bhat
        } else {
          pred_list[[i]] = Xtest[, gamma] %*% Bhat
        }
        Yhtest = as.numeric(scale(Yhtest, center = FALSE))
        pred_list[[i]] = as.numeric(scale(pred_list[[i]], center = FALSE))
        sigmamat =  diag(length(pred_list[[i]])); diag(sigmamat) = sigmahat
        W[i] = 1 / dmvnorm(Yhtest, pred_list[[i]], sigma = sigmamat)
      }
      W = W / max(W)
      pred_list = lapply(1:length(pred_list), function(i) W[i] * pnorm(pred_list[[i]]))
      pred_list = Reduce("+", pred_list) / sum(W)
      Ypred[[k]][folds[[k]] == m] = pred_list
    }
  }
  return(Ypred)
}
