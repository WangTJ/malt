library(MALT)

p=3 ;q=4;m=2
v1 = c(0, 1, 3, 5)/2
v2 = v1*1.5
partition = c(1, 2, 2, 2, 2, 2, 3, 3, 4, 3, 3, 4)
# coefficient grid
grid1 = matrix(v1[partition], p , q)
grid2 = matrix(v2[partition], p , q)

# patient distribution
p1 = matrix(rep(1/(p*q) , p*q) , p , q)
p1 = outer(c(0.5, 0.3, 0.2), c(0.4, 0.3, 0.2, 0.1))
p2 = outer(c(0.2, 0.3, 0.5), c(0.1, 0.15, 0.25, 0.5))
p2 = outer(c(0.5, 0.5, 0.), c(0.25, 0.25, 0.25, 0.25))


beta.real = as.param(c(as.vector(grid1), as.vector(grid2)), p, q, m)
beta0.real = as.param(as.vector(grid1), p, q)

n1 = 500; n2 = 600;
n = n1 + n2

data1 = dataset.gen(grid1, p1 , n1, cen.prob = 0.2, setseed = TRUE)
data2 = dataset.gen(grid2, p2 , n2, cen.prob = 0.2, setseed = TRUE)

data1$"G" = rep("sim1", n1)
data2$"G" = rep("sim2", n2)
data3$"G" = rep("sim3", n3)
data = rbind(data1, data2)

inibeta0 <- seq(0, 5, length.out=p*q)
inibeta <- rep(inibeta0, m)

beta.cox <- lasso.tree(data$G, data$T, data$N, data$time, data$status, 0, inibeta, maxiter = 100, eps = 1e-4)
beta.out <- lasso.tree(data$G, data$T, data$N, data$time, data$status, 4.6, inibeta, maxiter = 100, eps = 1e-3)

plot(beta.out)
plot(beta.cox)

lambda <- n*exp(seq(-7.5, -3.3, length.out = 20))
beta.seq <- lasso.tree(data$G, data$T, data$N, data$time, data$status, lambda, inibeta, maxiter = 100, eps = 1e-3)
beta.bic <- lasso.tree.bic(data$G, data$T, data$N, data$time, data$status, lambda, inibeta, maxiter = 100, eps = 1e-2)

plot(lambda, beta.bic$bics)
plot(lambda, similarity(beta.real, beta.bic$beta.seq)$similarity)
plot(as.param(beta.bic$beta.seq[1, ], p, q, m))
plot(as.param(beta.bic$beta.seq[20, ], p, q, m))
plot(beta.bic$beta.opt)
PlotCurve(beta.bic$beta.opt, data$G, data$T, data$N, data$time, data$status, ncol = 2)

similarity(beta.bic$beta.opt , beta.real)
