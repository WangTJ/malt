#functions for simulations and testing

#'Generate one dataset
#'
#'@param beta.mat beta's in a matrix form
#'@param n sample size
#'@param setseed boolen variable specifying whether fix seed
#'@param cen.prob probability of censoring
#'
#'@return data.frame of one dataset
#'
#'@export
dataset.gen <- function(beta.mat, p.mat, n, setseed = TRUE, cen.prob = 0){
  if (setseed) set.seed(123)
  p <- dim(beta.mat)

  TN <- sample(0:(p[1]*p[2]-1) ,n , replace = TRUE , prob = as.vector(p.mat))

  T <- TN %/% p[2] + 1
  N <- TN %% p[2] + 1


  y = rexp(n, rate = exp(beta.mat[cbind(T,N)]))
  eta <-  seq(from = min(beta.mat)/2 , to = max(beta.mat)*2 , length.out = 200)

  targetdelta <- rep(0, n)
  for (i in 1:length(eta))
  {
    # Exponential Cencoring
    C <- rexp(n, rate = exp(eta[i]))
    tilda_y <- pmin(C,y)
    delta <- as.numeric(y<=C)
    if (abs(1-mean(delta) - cen.prob)<abs(1-mean(targetdelta) - cen.prob)) {
      targety <- tilda_y; targetdelta <- delta;
    }
  }
  data.frame(T = T, N = N, time = targety, status = targetdelta)
}
