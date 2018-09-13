library(MALT)
library(parallel)

Scenario1 <- function(i) {
  print(i)
  set.seed(i)
  p=3 ;q=4;m=2
  v1 = c(0, 1, 3, 5)/2
  v2 = v1*1.5
  partition = c(1, 2, 2, 2, 2, 2, 3, 3, 4, 3, 3, 4)
  # coefficient grid
  grid1 = matrix(v1[partition], p , q)
  grid2 = matrix(v2[partition], p , q)

  # patient distribution
  #p1 = matrix(rep(1/(p*q) , p*q) , p , q)
  p1 = outer(c(0.5, 0.3, 0.2), c(0.4, 0.3, 0.2, 0.1))
  p2 = outer(c(0.1, 0.3, 0.6), c(0.1, 0.2, 0.3, 0.4))
  #p2 = outer(c(0.5, 0.5, 0.), c(0.25, 0.25, 0.25, 0.25))


  beta.real = as.param(c(as.vector(grid1), as.vector(grid2)), p, q, m)
  beta0.real = as.param(as.vector(grid1), p, q)

  n1 = 600; n2 = 600;
  n = n1 + n2

  data1 = dataset.gen(grid1, p1 , n1, cen.prob = 0.2, setseed = TRUE)
  data2 = dataset.gen(grid2, p2 , n2, cen.prob = 0.2, setseed = TRUE)

  data1$"G" = rep("trial1", n1)
  data2$"G" = rep("trial2", n2)

  data = rbind(data1, data2)

  inibeta0 <- seq(0, 5, length.out=p*q)
  inibeta <- rep(inibeta0, m)
  lambda <- n*exp(seq(-10, -2, length.out = 20))

  # MALT Result
  beta.meta.bic <- lasso.tree.bic(data$G, data$T, data$N, data$time, data$status, lambda, inibeta, maxiter = 100, eps = 1e-2)
  meta.result = similarity(beta.real, beta.meta.bic$beta.opt )

  # Single Result
  beta.trial1.bic <- lasso.tree.bic(data1$G, data1$T, data1$N, data1$time, data1$status, lambda, inibeta0, maxiter = 100, eps = 1e-2)
  trial1.result = similarity(beta0.real, beta.trial1.bic$beta.opt )

  beta.trial2.bic <- lasso.tree.bic(data2$G, data2$T, data2$N, data2$time, data2$status, lambda, inibeta0, maxiter = 100, eps = 1e-2)
  trial2.result = similarity(beta0.real, beta.trial2.bic$beta.opt )

  # Pooling Result
  data$"G" = rep("pooling" , n1+n2)
  beta.meta.bic <- lasso.tree.bic(data$G, data$T, data$N, data$time, data$status, lambda, inibeta0, maxiter = 100, eps = 1e-2)
  pooling.result = similarity(beta0.real, beta.meta.bic$beta.opt )

  return(list(meta.sse=as.numeric(meta.result$sse) , meta.sim=meta.result$similarity ,
              trail1.sse=as.numeric(trial1.result$sse), trial1.sim=trial1.result$similarity ,
              trial2.sse=as.numeric(trial2.result$sse), trial2.sim=trial2.result$similarity ,
              pooling.sse=as.numeric(pooling.result$sse), pooling.sim=pooling.result$similarity))
}

Result <- t(simplify2array(mclapply(1:8, Scenario1,mc.cores=8L)))

