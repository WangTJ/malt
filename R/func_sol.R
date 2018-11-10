#functions to solve the lasso tree problem


#'The core function for lasso-tree with meta-analysis
#'
#'Find the best estimate of parameter beta's, for tuning parameters lambda
#'
#'@param G label of dataset
#'@param T,N features one of the tumor
#'@param y time to failure
#'@param cen censor indicator
#'@param lambda tuning parameter of the fused group lasso penalty, could be a single value or a vector
#'@param inibeta vector of initial guess of beta's
#'@param trace boolen varable whether to show the process of calculation
#'@param maxiter max number of iteration
#'@param eps tolerance
#'@param w_select the method to select the adaptive weight w, "plainCox" or "preDefine"
#'@param w adaptive weight w
#'
#'@import survival
#'
#'@return the estimate of lambda in "param" class (see \code{\link{as.param}} for details)
#'
#'@examples
#'
#'@export
lasso.tree <- function(G, T, N, y, cen, lambda, inibeta = NULL, trace=TRUE, maxiter=200, eps=1e-4, w_select="plainCox", w=NULL) {
  n <- length(T)
  m <- length(unique(G))
  T <- as.factor(T)
  N <- as.factor(N)
  G <- as.factor(G)
  p <- length(unique(T))
  q <- length(unique(N))
  nb <- 2*p*q - p - q #number of boundaries
  G.levels <- levels(G)

  #create the design matrix x
  datap <- data.frame(y)
  datap$TN <- interaction(T, N)
  x <- model.matrix(~TN, datap)
  id <- which(apply(x, 1, sum) == 2)
  x[id, 1] <- 0
  colnames(x)[1] <- "TN1.1"

  #create a matrix g.index as the dataset indicator
  g.index = NULL
  for (o in levels(G)) g.index <- cbind(g.index, G == o)

  #the design matrix X to be multiplied by the beta_vec = (beta1, beta2, ...)'
  X <- matrix(0, n, p*q*m)
  for (o in 1:m) X[g.index[, o], (p*q*(o-1)+1):(p*q*o)] <- x[g.index[, o], ]

  ## calculate risk set
  indfail <- seq(1, n)[cen == 1]
  indset <- matrix(F, nrow=n, ncol=length(indfail))
  for (j in 1:length(indfail))
    indset[,j] <- (y >= y[indfail[j]]) & (G == G[indfail[j]]) #risk set is defined within each data set

  ## Create matrix of the partial ordering constraints for one dataset
  Apo0 <- matrix(0, nb, p*q)
  if(q > 1)
    for(i in 1:(p*(q - 1))) {
      Apo0[i, i] <- -1
      Apo0[i, i + p] <- 1
    }
  if(p > 1)
    for(i in 1:(q*(p - 1))) {
      j <- (i - 1)%/%(p-1) + i
      Apo0[i + p*(q - 1), j:(j + 1)] <- c(-1, 1)
    }

  #create matrix of partial ordering constraints for all datasets
  Apo <- Apo0
  if(m>1)
    for(o in 1:(m-1)) Apo <- Matrix::bdiag(Apo, Apo0)
  Apo <- as.matrix(Apo)

  ## Fitting the plain cox model to get w
  if (w_select == "plainCox") {
    suppressWarnings(model.cox <- coxph(Surv(y,cen)~X))
    #Apo %*% as.vector(model.cox$coefficients)
    w = NULL
    for (i in 1:(nb*m))
      w = c(w , Apo[i,which(Apo[i,]!=0)] %*% as.vector(model.cox$coefficients)[which(Apo[i,]!=0)])
    w[is.na(w)]=0
    w = pmax(w,0)
    w = sqrt(rowSums(matrix(w^2,nb,m))+0.0001)^(-1)
  }


  if (missing(inibeta)) inibeta <- rep(seq(0, 1, length.out = p*q), m)

  if (length(lambda) == 1){
    beta.out <- lasso_tree_single(X, indset, Apo, nb, m, cen, lambda, w, inibeta, maxiter, eps, trace)
  }
  else{
    beta.out <- lasso_tree_multi(X, indset, Apo, nb, m, cen, lambda, w, inibeta, maxiter, eps, trace)
  }
  beta.out <- as.param(beta.out, p, q, m, G.levels, lambda)
  return(beta.out)
}


#'Function to calculate the negative log partial likelihood
#'
#'@param G label of dataset
#'@param T,N features one of the tumor
#'@param y time to failure
#'@param cen censor indicator
#'@param beta vector of parameter beta's
#'
#'@return the negative log partial likelihood
#'
#'@export
logPL <- function(G, T, N, y, cen, beta)
{
  n <- length(T)
  m <- length(unique(G))
  T <- as.factor(T)
  N <- as.factor(N)
  G <- as.factor(G)
  p <- length(unique(T))
  q <- length(unique(N))

  #create the design matrix x
  datap <- data.frame(y)
  datap$TN <- interaction(T, N)
  x <- model.matrix(~TN, datap)
  id <- which(apply(x, 1, sum) == 2)
  x[id, 1] <- 0
  colnames(x)[1] <- "TN1.1"

  #create a matrix g.index as the dataset indicator
  g.index = NULL
  for (o in levels(G)) g.index <- cbind(g.index, G == o)

  #the design matrix X to be multiplied by the beta_vec = (beta1, beta2, ...)'
  X <- matrix(0, n, p*q*m)
  for (o in 1:m) X[g.index[, o], (p*q*(o-1)+1):(p*q*o)] <- x[g.index[, o], ]

  ## calculate risk set
  indfail <- seq(1, n)[cen == 1]
  indset <- matrix(F, nrow=n, ncol=length(indfail))
  for (j in 1:length(indfail))
    indset[,j] <- (y >= y[indfail[j]]) & (G == G[indfail[j]]) #risk set is defined within each data set

  return(logPL_C(X, indset, cen, beta))
}

#'Function to select the best beta according to BIC
#'
#'For a vector of tuning parameters lambda, calculate all the beta and select the best one according to
#'  BIC criterion
#'
#'@param G label of dataset
#'@param T feature one of the tumor
#'@param N feature two of the tumor
#'@param y time to failure
#'@param cen censor indicator
#'@param lambda tuning parameter of the fused group lasso penalty, a vector
#'@param inibeta vector of initial values of beta's
#'@param trace boolen varable whether to show the process of calculation
#'@param maxiter max number of iteration
#'@param eps tolerance
#'@param w_select the method to select the adaptive weight w, "plainCox" or "preDefine"
#'@param w adaptive weight w
#'
#'@return a list including: the best solution of beta according to BIC, a vector of BIC values corresponding
#'  to the input lambda, and a "param" beta.seq including all the beta's
#'
#'@examples
#'
#'@export
lasso.tree.bic <- function(G, T, N, y, cen, lambda, inibeta = NULL, trace=TRUE, maxiter=200, eps=1e-4, w_select="plainCox", w=NULL){
  nlambda <- length(lambda)
  p <- length(unique(T))
  q <- length(unique(N))
  m <- length(unique(G))
  G.levels <- levels(as.factor(G))

  beta.seq <- lasso.tree(G, T, N, y, cen, lambda, inibeta, trace, maxiter, eps, w_select, w)
  bics <- bic(beta.seq, G, T, N, y, cen)

  ind.opt <- which.min(bics)
  beta.opt <- as.param(beta.seq[ind.opt, ], p, q, m, G.levels, lambda[ind.opt])

  return(list(beta.opt = beta.opt, bics = bics, beta.seq = beta.seq))
}

