setwd(this.path::here())
# Case 1: all entries of beta are positive

#######################################################
# Case 1: positive betas, transformed region
######################################################

library(mig)
set.seed(12345)
beta <- rexp(2)
p1 <- q <- rexp(2)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
pts <- rbind(p1, p2, p3)
pdf("../_manuscript/figures/pmig_rotation.pdf", width = 6, height = 9)
par(mfrow = c(3,2), mar = c(4,4,1,1), bty = "l")
yran <- range(pts[,2]) + c(-0.1, 0.2)
xran <- range(pts[,1]) + c(-0.1, 0.2)
plot(NA, xlab = expression(X[1]),
     ylab = expression(X[2]),
     ylim = yran,
     xlim = xran,
     yaxs = "i",
     xaxs = "i",
     panel.first = abline(a=0, b = -beta[1]/beta[2]))
polygon(c(-yran[1]*beta[2]/beta[1], xran[1], xran[1]),
        c(yran[1], -beta[1]/beta[2]*xran[1], yran[1]),
        col = "grey90")
polygon(pts, col = 4)

set.seed(12345)
beta <- rexp(2)
p1 <- q <- rexp(2)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
pts <- rbind(p1, p2, p3)

Q2 <- c(-beta[2], beta[1])/sqrt(sum(beta^2))
Q <- rbind(beta, Q2)
tcrossprod(Q)
trpts <- pts %*% t(Q)
plot(NA, xlab = "R",
     ylab = "Z",
     xaxs = "i",
     ylim = range(trpts[,2]),
     xlim = range(trpts[,1]))
# polygon(x = c(-1, 0, 0, -1),
        # y = rep(c(range(trpts[,2])) + c(-1,1), each = 2), col = "grey90")
polygon(trpts, col = 4)
v1 <- q
v2 <- c(q[1], - beta[1]*q[1]/beta[2])
v3 <- c(-beta[2]*q[2]/beta[1], q[2])

zmin = Q2 %*% v2
zmax = Q2 %*% v3





#######################################################
# Case 2: betas of opposite signs, beta[1] negative
######################################################


# Case 2: betas of opposite signs, beta[1] negative
set.seed(1234)
beta <- c(-rexp(1)/4 -1 , rexp(1)/4 + 1)
xi <- beta + runif(2)*sign(beta)
Omega <- diag(2)

p1 <- q <- beta + runif(2, 0, 0.1)*sign(beta)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
p4 <- c(-Inf, -Inf)
pts <- rbind(p1, p2, p3, p4)
pts <- pts[rowSums(scale(pts, center = q, scale = FALSE) <= 0) == 2L,]
samp <- rmig(n = 2500, xi = xi, Omega = Omega, beta = beta)
# par(mfrow = c(1,2))
plot(NA, xlab = expression(X[1]),
     ylab = expression(X[2]),
     ylim = c(-10,5),
     xlim = c(-10,10),
     xaxs = "i",
     yaxs = "i",
     panel.first = abline(a = 0, b = -beta[1]/beta[2]))
# points(samp, pch = 20, col = scales::alpha("black", alpha = 0.1))

minpts <- apply(samp, 2, min) - 100
modifpts <- rbind(pts[1:2,],
                  # The follow two points are substitutes for -(Inf, Inf)
                  c(-beta[2]*minpts[2]/beta[1], minpts[2]),
                  minpts,
                  c(minpts[1] - 2, q[2]),
                  c(-15, q[2]))
polygon(x = c(beta[2]*10/beta[1], 10, 10, -5*beta[2]/beta[1]),
        y = c(-10, -10, 5, 5),
        col = "grey90")
polygon(modifpts, col = 4)


set.seed(1234)
beta <- c(-rexp(1)/4 -1 , rexp(1)/4 + 1)
xi <- beta + runif(2)*sign(beta)
Omega <- diag(2)

p1 <- q <- beta + runif(2, 0, 0.1)*sign(beta)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
p4 <- c(-Inf, -Inf)
pts <- rbind(p1, p2, p3, p4)
pts <- pts[rowSums(scale(pts, center = q, scale = FALSE) <= 0) == 2L,]
samp <- mig::rmig(n = 2500, xi = xi, Omega = Omega, beta = beta)

minpts <- apply(samp, 2, min) - 100
modifpts <- rbind(pts[1:2,],
                  # The follow two points are substitutes for -(Inf, Inf)
                  c(-beta[2]*minpts[2]/beta[1], minpts[2]),
                  minpts,
                  c(minpts[1] - 2, q[2]),
                  c(-15, q[2]))
Q2 <- c(-beta[2], beta[1])/sqrt(sum(beta^2))
Q <- rbind(beta, Q2)
tcrossprod(Q)
trpts <- pts %*% t(Q)
 trpts <- modifpts %*% t(Q)
plot(NA, xlab = "R",
     ylab = "Z",
     ylim = c(-1,6),
     xlim = c(0,5),
     xaxs = "i")
polygon(trpts, col = 4)
points(beta %*% q, Q2 %*% q)
tq <- Q %*% q
v1 <- c(q[1], - beta[1]*q[1]/beta[2])
v2 <- c(-beta[2]*q[2]/beta[1], q[2])
tv1 <- Q %*% v1
tv2 <- Q %*% v2
z_r <- function(r){
        pmax(tv1[2] + Q2[2]/beta[2]*r, tv2[2] + Q2[1]/beta[1]*r)}



#########################################################
## Case 3, both betas negative
#########################################################


set.seed(1234)
beta <- (-rexp(2) - 2)/4
xi <- beta + runif(2)*sign(beta)
Omega <- diag(2)

p1 <- q <- beta + runif(2, 0, 0.1)*sign(beta)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
p4 <- c(-Inf, -Inf)
pts <- rbind(p1, p2, p3, p4)
pts <- pts[rowSums(scale(pts, center = q, scale = FALSE) <= 0) == 2L,]
samp <- rmig(n = 2500, xi = xi, Omega = Omega, beta = beta)
# par(mfrow = c(1,2))
plot(NA, xlab = expression(X[1]),
     ylab = expression(X[2]),
     ylim = c(-10,5),
     xlim = c(-5,5),
     xaxs = "i",
     yaxs = "i",
     panel.first = abline(a=0, b = -beta[1]/beta[2]))

minpts <- apply(samp, 2, min) - 100
modifpts <- rbind(pts[1,],
                  c(q[1], minpts[2]),
                  minpts,
                  c(minpts[1], q[2]))
polygon(x = c(5, 5,  -beta[2]/beta[1]*5),
        y = c(-beta[1]/beta[2]*5, 5, 5),
        col = "grey90")
polygon(modifpts, col = 4)

set.seed(1234)
beta <- (-rexp(2) - 2)/4
xi <- beta + runif(2)*sign(beta)
Omega <- diag(2)

p1 <- q <- beta + runif(2, 0, 0.1)*sign(beta)
p2 <- c(q[1], -beta[1]*q[1]/beta[2])
p3 <- c(-beta[2]*q[2]/beta[1], q[2])
p4 <- c(-Inf, -Inf)
pts <- rbind(p1, p2, p3, p4)
pts <- pts[rowSums(scale(pts, center = q, scale = FALSE) <= 0) == 2L,]
samp <- rmig(n = 2500, xi = xi, Omega = Omega, beta = beta)

minpts <- apply(samp, 2, min) - 100
modifpts <- rbind(pts[1,],
                  c(q[1], minpts[2]),
                  minpts,
                  c(minpts[1], q[2]))

Q2 <- c(-beta[2], beta[1])/sqrt(sum(beta^2))
Q <- rbind(beta, Q2)
tcrossprod(Q)
trpts <- modifpts %*% t(Q)
plot(NA, xlab = "R",
     ylab = "Z",
     ylim = c(-1,5),
     xlim = c(0,7),
     xaxs = "i")
polygon(trpts, col = 4)

dev.off()
