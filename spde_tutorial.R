# In the literature there is some ideas to fit GF by an approximation of the GF to any
# GMRF. The very good alternative found is shuch one that found a explicit link between
# the an stochastic partial differential equation (SPDE) and the Mat´ern Gaussian random
# fields, [Lindgren et al., 2011]. The solution of the that SPDE thought the Finite Element
# Method (FEM) provides a explicit link between GF to GMRF.
#
# Code from Chapter 8 of http://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/spde-tutorial.pdf
#
# LGCP example
require(spatstat)
win <- owin(c(0,3), c(0,3))

spatstat.options(npixel=300)
beta0 <- 3
# expected number of points
print(exp(beta0) * diff(range(win$x)) * diff(range(win$y)))

sigma2x <- 0.2; kappa <- 2

# do simulation
set.seed(0)
lg.s <- rLGCP('matern', beta0, var=sigma2x, scale=1/kappa, nu=1, win=win)

# point pattern locations
print(n <- nrow(xy <- cbind(lg.s$x, lg.s$y)[,2:1]))

Lam <- attr(lg.s, 'Lambda')
print(summary(as.vector(rf.s <- log(Lam$v))))

par(mfrow=c(1,1))
require(fields)
image.plot(list(x=Lam$yrow, y=Lam$xcol, z=rf.s), main='log-Lambda', asp=1)
points(xy, pch=19)

# create a mesh for inference
loc.d <- 3*t(matrix(c(0,0,1,0,1,1,0,1,0,0), 2))
(nv <- (mesh <- inla.mesh.2d(loc.d=loc.d, off=.2, max.e=.5, cut=.1))$n)

par(mar=c(0,0,0,0))
plot(mesh, asp=1, main='')
points(xy, col=4, pch=19); lines(loc.d, col=3)

# define SPDE model
spde <- inla.spde2.matern(mesh=mesh, alpha=2)

# but mesh has nodes out of the domain, so the mesh area is too big:
print(sum(diag(spde$param.inla$M0)))

# build mesh without outer extension
require(deldir)
dd <- deldir(mesh$loc[,1], mesh$loc[,2])
tiles <- tile.list(dd)
require(gpclib)
area.poly(pl.study <- as(loc.d, 'gpc.poly'))
sum(w <- sapply(tiles, function(p) area.poly(intersect(as(cbind(p$x, p$y), 'gpc.poly'), pl.study))))
print(table(w>0))

par(mar=c(2,2,1,1), mgp=2:0)
plot(mesh$loc, asp=1, col=(w==0)+1, pch=19, xlab='', ylab='')
for (i in 1:length(tiles))
  lines(c(tiles[[i]]$x, tiles[[i]]$x[1]), c(tiles[[i]]$y, tiles[[i]]$y[1]))
lines(loc.d, col=3)

# the data
y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n))

# projector matrix for observed location points and integration points
lmat <- inla.spde.make.A(mesh, xy)
imat <- Diagonal(nv, rep(1, nv))
A.pp <- rBind(imat, lmat)
stk.pp <- inla.stack(data=list(y=y.pp, e=e.pp),
                     A=list(1,A.pp), tag='pp',
                     effects=list(list(b0=rep(1,nv+n)), list(i=1:nv)))

# run inla
pp.res <- inla(y ~ 0 + b0 + f(i, model=spde),
               family='poisson', data=inla.stack.data(stk.pp),
               control.predictor=list(A=inla.stack.A(stk.pp)),
               E=inla.stack.data(stk.pp)$e)
pp.rf <- inla.spde2.result(pp.res, 'i', spde)

# plot posteriors
par(mfrow=c(2,2), mar=c(3,3,1,0.3), mgp=c(2,1,0))
plot(pp.res$marginals.fix[[1]], type='l',
     xlab=expression(beta[0]), ylab='Density')
abline(v=beta0, col=2)
plot(pp.rf$marginals.variance.nominal[[1]], type='l',
     xlab=expression(sigma^2), ylab='Density')
abline(v=sigma2x, col=2)
plot(pp.rf$marginals.kappa[[1]], type='l',
     xlab=expression(kappa), ylab='Density')
abline(v=kappa, col=2)
plot(pp.rf$marginals.range.nominal[[1]], type='l',
     xlab='Nominal range', ylab='Density')
     abline(v=sqrt(8*1)/kappa, col=2)