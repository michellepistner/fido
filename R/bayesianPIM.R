GG <- function(obj, D){
  if(is.null(nrow(obj)) | nrow(obj) != ncol(obj)){
    F = cbind(diag(D-1), rep(-1,D-1))
    H = matrix(1,1,D)
    
    obj.par = F%*% obj
    obj.perp = H %*% obj
    return(list(obj.par = obj.par, obj.perp = obj.perp))
  } else{
    F = cbind(diag(D-1), rep(-1,D-1))
    H = matrix(1,1,D)
    
    obj.par = F%*% obj %*% t(F)
    obj.perp = H %*% obj %*% t(H)
    obj.plus = H %*% obj %*% t(F)
    return(list(obj.par = obj.par, obj.perp = obj.perp, obj.plus = obj.plus))
  }
}

t.sampler <- function(nu.star, M.star, Xi.star, V.star){
  Sigma = rinvwishart(nu.star, Xi.star)
  C = t(chol(Sigma))
  X = rmatnorm(1,M.star,diag(nrow(M.star)), V.star)
  Y= C %*% X
  
  return(Y)
}


bayesPIM_scaleModel <- function(fitc){
  n_samples <- dim(fitc$Samples)[3]
  N <- dim(fitc$Samples)[2]
  D <- dim(fitc$Samples)[1] + 1 ##ALR coords so add 1
  tau = matrix(NA, nrow = N, ncol = n_samples)
  
  lambda.par = alrInv_array(fitc$Samples, coords = 1)
  FF = cbind(diag(D-1), -1)
  lambda.par_tmp = array(NA, dim = c(dim(lambda.par)[1] - 1, dim(lambda.par)[2:3]))
  for(i in 1:n_samples){
    lambda.par_tmp[,,i] = FF %*% lambda.par[,,i]
  }
  
  trans_priors <- transformedPIM_priors(Theta, Xi)
  Xi.t = trans_priors$Xi
  Theta.t = trans_priors$Theta
  
  Theta.trans = GG(Theta.t, nrow(Xi.t))
  Xi.trans = GG(Xi.t, nrow(Xi.t))
  
  FF = cbind(diag(D-1), rep(-1,D-1))
  for(i in 1:n_samples){
    nu.star = upsilon + D
    M.star = Theta.trans$obj.perp %*% X +  Xi.trans$obj.plus %*% solve(Xi.trans$obj.par) %*% ((lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X))
    V.star = (diag(N) + t(X) %*% Gamma %*% X) %*% (diag(N) + solve((diag(N) + t(X) %*% Gamma %*% X)) %*% (t(lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X) %*% solve(Xi.trans$obj.par) %*% ((lambda.par_tmp[,,i] - Theta.trans$obj.par %*% X))))
    Xi.star = Xi.trans$obj.perp - (Xi.trans$obj.plus) %*% solve(Xi.trans$obj.par) %*% t(Xi.trans$obj.plus)
    tau[,i] = t.sampler(nu.star, M.star, Xi.star, V.star)
  }
  
  return(list(tau=tau))
}

supp_wTotals <- function(fitc, lambda.total){
  n_samples <- dim(fitc$Samples)[3]
  W.par <- alrInv_array(fitc$Samples, coords = 1)
  lambda.par <- log(W.par)
  lambda <- array(NA, dim = dim(lambda.par)) #Always needed
  for(i in 1:n_samples){
    lambda[,,i]  <- sweep(lambda.par[,,i], MARGIN=2, lambda.total$tau[,i], `+`)
  }
  ##Transform samples back
  return(lambda)
}

transformedPIM_priors <- function(theta.old, Xi.old){
  
  theta.new <- alrInv_array(theta.old, D, 1)
  
  ##Last piece: how to transform Xi
  ##Xi.old --> G Xi G^T
  ##Xi.old^-1 -->
  Xi.new <- Xi.old
  return(list(Theta = theta.new, Xi = Xi.new))
}