#functions for visualization

#'Plot the beta matrix as tiles
#'
#'This function plot the values of beta's in a tile plot: use T and N as axes, and color to represent the values
#'  of beta's.
#'
#'@param beta.mat the value of beta's in a matrix form
#'
#'@import ggplot2
#'@import ggrepel
#'
#'@export
tile.plot <- function(beta.mat)
{
  nT <- dim(beta.mat)[1]
  nN <- dim(beta.mat)[2]
  beta <- as.vector(beta.mat)
  T <- rep(seq(1, nT), nN)
  N <- NULL
  for (i in 1:nN) N <- c(N, rep(i, nT))
  betaframe <- data.frame(T = T, N = N, beta = beta)
  ggplot(betaframe, aes(x = N, y = T)) + geom_raster(aes(fill = beta)) + coord_fixed(ratio = 1)
}

#'Plot function for class "param"
#'
#'This function plot the coefficients as an object "param" (see \code{\link{as.param}} for details) on a 2-way
#'  T/N table, using the darkness of the color to represent the values of coefficients.
#'
#'@param beta beta is of class "param"
#'@param lambda.index when the param beta contains multiple values of the tuning parameter, the user need to
#'  specify the index of the tuning parameter to be plotted
#'
#'@import ggplot2
#'
#'@export
plot.param <- function(beta, lambda.index = NULL)
{
  lambda <- attr(beta, "lambda")
  if(length(lambda)>1 & is.null(lambda.index)){
    cat("There are", length(lambda), "tuning parameters:", lambda, "\n", "Please specify a index.\n")
  }else{
    p <- attr(beta, "p")
    q <- attr(beta, "q")
    m <- attr(beta, "m")
    G.levels <- attr(beta, "Glevels")
    if(is.null(G.levels)) G.levels <- 1:m

    T <- rep(rep(seq(1, p), q), m)

    N <- NULL
    for (i in 1:q) N <- c(N, rep(i, p))
    N <- rep(N, m)

    G <- NULL
    for (i in 1:m) G <- c(G, rep(G.levels[i], p*q))

    if(length(lambda) <= 1){
      betaframe <- data.frame(G = G, T = T, N = N, beta = as.vector(beta))
    }else{
      betaframe <- data.frame(G = G, T = T, N = N, beta = as.vector(beta[lambda.index, ]))
    }

    ggplot(betaframe, aes(x = N, y = T)) + geom_raster(aes(fill = beta)) + coord_fixed(ratio = 1) + facet_wrap(~G)
  }

}

#'Constructor of an object "param", which is used for better storage, printing and visualization of the coefficients.
#'
#'This function constructs a "param" beta from a vector or a matrix of values of the coefficients. A param object
#'  can have multiple values of tuning parameter lambda.
#'
#'@param beta.data the values of the coefficients. For a single tuning parameter, it is a vector; for multiple
#'  tuning parameters, it is a matrix with each row associated with one tuning parameter
#'@param p number of levels of T
#'@param q number of levels of N
#'@param m number of data sets (levels of G)
#'@param G.levels levels of group indicator
#'@param lambda the tuning parameter associated with this beta, could be a scalar or a vector
#'
#'@export
as.param <- function(beta.vec, p, q, m = 1, G.levels = NULL, lambda = NULL)
{
  beta.param <- beta.vec
  #force the first beta of each data set to be 0
  if(length(lambda) <= 1){ #if lambda = NULL, then length(lambda) = 0
    for(o in 1:m) beta.param[(p*q*(o-1)+1):(o*p*q)] <- beta.param[(p*q*(o-1)+1):(o*p*q)] - beta.param[p*q*(o-1)+1]
  }else{
    for(o in 1:m) beta.param[ , (p*q*(o-1)+1):(o*p*q)] <- beta.param[ , (p*q*(o-1)+1):(o*p*q)] - beta.param[ , p*q*(o-1)+1]
  }
  attr(beta.param, "p") <- p
  attr(beta.param, "q") <- q
  attr(beta.param, "m") <- m
  attr(beta.param, "Glevels") <- G.levels
  attr(beta.param, "lambda") <- lambda
  class(beta.param) <- "param"
  return(beta.param)
}

#'Summary of the class "param"
#'
#'This function prints out a "param" object (see \code{\link{as.param}} for details) in a decent way. It only
#'  works for a param with a single lambda value.
#'
#'@param beta an object in class "param"
#'@param lambda.index when the param beta contains multiple values of the tuning parameter, the user need to
#'  specify the index of the tuning parameter to be summaried
#'
#'@return a matrix describe he combined staging rules
#'@export
summary.param <- function(beta, lambda.index = NULL){
  lambda <- attr(beta, "lambda")
  if(length(lambda)>1 & is.null(lambda.index)){
    cat("There are", length(lambda), "tuning parameters:", lambda, "\n", "Please specify a index.\n")
  }else{
    p <- attr(beta, "p")
    q <- attr(beta, "q")
    m <- attr(beta, "m")
    G.levels <- attr(beta, "Glevels")
    if(is.null(G.levels)) G.levels <- 1:m
    Tnames <- paste("T", 1:p, sep = "")
    Nnames <- paste("N", 0:(q-1), sep = "")
    cat("lambda = ", lambda, "\n")
    for(o in 1:m){
      cat("Group ", G.levels[o], "\n")
      if(length(lambda)<=1) beta0 <- matrix(beta[((o-1)*p*q+1):(o*p*q)], p, q)
      else beta0 <- matrix(beta[lambda.index, ((o-1)*p*q+1):(o*p*q)], p, q)
      rownames(beta0) <- Tnames
      colnames(beta0) <- Nnames
      print(beta0)
    }

    nStage <- 1
    if (length(lambda)<=1){
      beta.stage <- (beta[1:(p*q)])^2
      for(o in 2:m){
        beta.stage <- beta.stage + (beta[((o-1)*p*q+1):(o*p*q)])^2
      }
      beta.stage <- sqrt(beta.stage)
    }
    else
    {
      beta.stage <- (beta[lambda.index,1:(p*q)])^2
      for(o in 2:m){
        beta.stage <- beta.stage + (beta[lambda.index,((o-1)*p*q+1):(o*p*q)])^2
      }
      beta.stage <- sqrt(beta.stage)
    }
    nStage <- length(unique(beta.stage))
    beta.stage <- as.factor(beta.stage)
    levels(beta.stage) = 1:nStage
    beta.stage = as.integer(beta.stage)
    beta.stage <- matrix(beta.stage, p, q)
    rownames(beta.stage) <- Tnames
    colnames(beta.stage) <- Nnames
    cat("Combined staging rules \n")
    print(beta.stage)
    return(beta.stage)
  }
}

#'Function to calculate the degrees of freedom for a "param" beta
#'
#'@param beta beta of class "param"
#'
#'@return degrees of freedom of each data set, a vector when beta includes one lambda, a matrix when beta
#'  includes multiple lambda
#'
#'@export
DoF.param <- function(beta){
  p <- attr(beta, "p")
  q <- attr(beta, "q")
  m <- attr(beta, "m")
  nlambda <- length(attr(beta, "lambda"))
  if(nlambda <= 1){
    dof <- NULL
    for (o in 1:m) dof = c(dof, 2+length(unique(beta[(p*q*(o-1)+1):(p*q*o)])))
  }else{
    dof <- NULL
    for (j in 1:nlambda){
      d <- NULL
      for (o in 1:m) d = c(d, 2+length(unique(beta[j, (p*q*(o-1)+1):(p*q*o)])))
      dof <- rbind(dof, d)
    }
  }
  return(dof)
}

#'Calculate the BIC value for a centain beta
#'
#'This function calculate the BIC value from beta and the data. Here we use the number of unique beta's as the degree of freedom, and the
#'  number of patients of all the data sets as the sample size.
#'
#'@param beta an object "param"
#'@param G group indicator vector
#'@param T a vector of cancer feature T
#'@param N a vector of cancer feature N
#'@param y a vector of time to event
#'@param cen censor indicator vector
#'
#'@return the BIC value. If beta contains single lambda value, the function returns a value of bic; if the beta contains
#'  vector of lambda's, the function returns a vector of bic values.
#'
#'@export
bic <- function(beta, G, T, N, y, cen)
{
  n <- length(T)
  m <- length(unique(G))
  T <- as.factor(T)
  N <- as.factor(N)
  G <- as.factor(G)
  p <- length(unique(T))
  q <- length(unique(N))
  nlambda <- length(attr(beta, "lambda"))
  Gsize <- as.vector(table(G))

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
  X <- matrix(0, n, p*q*(m))
  for (o in 1:m) X[g.index[, o], (p*q*(o-1)+1):(p*q*o)] <- x[g.index[, o], ]

  ## calculate risk set
  indfail <- seq(1, n)[cen == 1]
  indset <- matrix(F, nrow=n, ncol=length(indfail))
  for (j in 1:length(indfail))
    indset[,j] <- (y >= y[indfail[j]]) & (G == G[indfail[j]])

  if(nlambda <= 1){
    ans <- 2*logPL_C(X, indset, cen, beta) + as.numeric(crossprod(DoF.param(beta), log(Gsize)))
  }else{
    ans <- NULL
    for (j in 1:nlambda){
      b <- as.param(beta[j, ], p, q, m)
      ans <- c(ans, 2*logPL_C(X, indset, cen, b))
    }
    ans <- ans + as.vector(DoF.param(beta) %*% log(Gsize))
  }
  return(ans)
}

#'Function to compare how close two "param" beta's are
#'
#'@param beta0 a "param" beta to be compared with, including one lambda value
#'@param beta a "param" beta that can have one or multiple lambda values
#'
#'@return a list containing: (1) the similarity between the two beta's, defined by the
#'  proportion of the correct boundaries; (2) sse, sum of squre error of all the differences
#'  on the boundary, with repect to beta0; (3) sim.vec, a vector indicating which boundaries
#'  are the same.
#'@export
similarity <- function(beta0, beta){
  p <- attr(beta, "p")
  q <- attr(beta, "q")
  m <- attr(beta, "m")
  nlambda <- length(attr(beta, "lambda"))
  nb <- 2*p*q - p - q #number of boundaries

  ## Create matrix of the partial ordering constraints for one dataset
  Apo0 <- matrix(0, nb, p*q)
  if(q>1)
    for(i in 1:(p*(q - 1))) {
      Apo0[i, i] <- -1
      Apo0[i, i + p] <- 1
    }
  if(p>1)
    for(i in 1:(q*(p - 1))) {
      j <- (i - 1)%/%(p-1) + i
      Apo0[i + p*(q - 1), j:(j + 1)] <- c(-1, 1)
    }
  #create matrix of partial ordering constraints for all datasets
  Apo <- Apo0
  if(m>1)
    for(o in 1:(m-1)) Apo <- Matrix::bdiag(Apo, Apo0)
  Apo <- as.matrix(Apo)

  bound0.val <- Apo %*% beta0
  bound0.ind <- (bound0.val) == 0

  if(nlambda <= 1){
    bound.val <- Apo %*% beta
    bound.ind <- (bound.val) == 0
    sim.vec <- as.vector(bound.ind == bound0.ind)
    sim <- sum(sim.vec)/(nb*m)
    sse <- crossprod(bound.val - bound0.val)
  }else{
    sim <- NULL
    sse <- NULL
    sim.vec <- NULL
    for(j in 1:nlambda){
      bound.val <- Apo %*% beta[j, ]
      bound.ind <- (bound.val) == 0
      sim.v <- as.vector(bound.ind == bound0.ind)
      sim <- c(sim, sum(sim.v))
      sim.vec <- rbind(sim.vec, sim.v)
      sse <- c(sse, crossprod(bound.val - bound0.val))
    }
    sim <- sim/(nb*m)
  }
  return(list(similarity = sim, sse = sse, sim.vec = sim.vec))
}

#'Function to plot the survival curves
#'
#'@param beta the estimated coefficients, in class "param"
#'@param G a vector of group labels
#'@param T,N cancer characteristic in vectors or factors
#'@param y time to events, in vector
#'@param cen censor indicator, in vector
#'@param ncol number of columns when plotting curves for each trial
#'@param legend.pos position of the legend, "bottomleft" or "topright"
#'
#'@import survival
#'
#'@export
PlotCurve <- function(beta, G, T, N, y, cen, ncol = 2, legend.pos = "bottomleft")
{
  p <- attr(beta, "p")
  q <- attr(beta, "q")
  m <- attr(beta, "m")
  G <- as.factor(G)
  Glevels <- levels(G)
  if(m != length(Glevels)) cat("Number of groups in beta and in G do not match!")

  nStage <- 1
  beta.stage <- (beta[1:(p*q)])^2

  for(o in 2:m){
    beta.stage <- beta.stage + (beta[((o-1)*p*q+1):(o*p*q)])^2
  }
  beta.stage <- sqrt(beta.stage)
  nStage <- length(unique(beta.stage))
  beta.stage <- as.factor(beta.stage)
  levels(beta.stage) = 1:nStage
  beta.stage = as.integer(beta.stage)
  beta.stage <- matrix(beta.stage, p, q)

  T.ind <- as.factor(T) #the indices of T starting from 1
  levels(T.ind) <- 1:(nlevels(T.ind))
  N.ind <- as.factor(N)
  levels(N.ind) <- 1:(nlevels(N.ind))
  stage <- beta.stage[cbind(T.ind, N.ind)]
  stage <- ordered(stage)
  df <- data.frame(G = G, time = y, cen = cen, stage = stage)

  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(ceiling(nlevels(G)/ncol), ncol))

  colorList = palette()
  for(trial in Glevels){
    df.sub <- df[df$G == trial, ]
    fit <- survfit(Surv(time, cen) ~ stage, data = df.sub)
    par(xpd = TRUE)
    plot(fit, mark.time = FALSE, lty = 1, col = colorList[sort(unique(df.sub$stage))], xlab = "Time",
         ylab = "Survival")
    #legend(max(df.sub$time), 1, levels(stage), col = 1:nlevels(stage), lty = 1, cex = 0.5)
    legend(legend.pos, levels(stage), col = colorList[1:nlevels(stage)], lty = 1, cex = 0.6, bty = "n")
    title(trial)
  }
  par(oldpar)
}
