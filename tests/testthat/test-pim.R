test_that("bayes PIM model works with no error!", {
  set.seed(1)
  sim <- pibble_sim()
  expect_error(expect_error(pibble(sim$Y, sim$X, isPIM = TRUE))) # expect no error!
})
# 
test_that("t sampler gives correct mean and covariance",{
  set.seed(1)
  n_samples <- 1000
  p <- 3
  m <- 2
  M <- matrix(1, nrow = 3, ncol = m)
  Sigma <- diag(p)
  Omega <- diag(m)
  t.star <- array(NA, c(dim(M), n_samples))
  for(i in 1:n_samples){
    t.star[,,i] = fido:::t.sampler(5,M,Sigma,Omega)
  }
  expect_true(sum(abs(apply(t.star,c(1,2),mean) - M )< 0.1) == (nrow(M)*ncol(M)))
  
  ##Testing covariance
  zz = rep(NA,n_samples)
  yy = rep(NA, n_samples)
  for(i in 1:n_samples){
    zz[i] = t.star[1,1,i]
    yy[i] = t.star[1,2,i]
  }
  ther.cov = (1/(9-2))* kronecker(Sigma, Omega)
  expect_true(abs(var(zz) - ther.cov[1,1]) < .1)
  expect_true(abs(cov(zz,yy) - ther.cov[2,1]) < .1)
})

