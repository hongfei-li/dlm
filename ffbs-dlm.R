################################################################################
### FFBS from TSA4
################################################################################

## packages
require(astsa)
require(truncnorm)

###-- Notation --###
## y(t) = x(t) + v(t); v(t) ~ iid N(0,V)
## x(t) = x(t-1) + w(t); w(t) ~ iid N(0,W)
## priors: x(0) ~ N(m0,C0); V ~ IG(a,b); W ~ IG(c,d)
## FFBS: x(t|t) ~ N(m,C); x(t|n) ~ N(mm,CC); x(t|t+1) ~ N(a,R) ##--

q025 <- function(x){quantile(x, 0.025)}
q975 <- function(x){quantile(x, 0.975)}
## FFBS
ffbs <- function(y, V, W, Phi, m0, C0){
    n  = length(y)
    a  = rep(0, n)
    R  = rep(0, n)
    m  = rep(0, n)
    C  = rep(0, n)
    B  = rep(0, n-1)
    H  = rep(0, n-1)
    mm = rep(0, n)
    CC = rep(0, n)
    x  = rep(0, n)
    llike = 0.0
    for (t in 1:n){
        ## forward filtering
        if(t==1){
            a[1] = Phi * m0 ## x_{t|t-1}
            R[1] = Phi * C0 * Phi + W ## P_{t|t-1}
        } else {
            a[t] = Phi * m[t-1]
            R[t] = Phi * C[t-1] * Phi + W
        }
        f = a[t] ## y_{t|t-1}
        Q = R[t] + V ## G_{t|t-1}
        A = R[t] / Q ## K_t: kalman gain
        m[t] = a[t] + A * (y[t] - f) ## x_{t|t}
        C[t] = R[t] - Q * A**2 ## P_{t|t}
        B[t-1] = Phi * C[t-1] / R[t] ## J_{t-1}
        H[t-1] = C[t-1] - B[t-1] * R[t] * B[t-1] ## G_{t-1|t-1}
        llike = llike + dnorm(y[t], f, sqrt(Q), log = TRUE)
    }
    mm[n] = m[n] ## sampled state mean
    CC[n] = C[n] ## sampled state variance
    x[n]  = rnorm(1, m[n], sqrt(C[n]))
    for (t in (n-1):1){
        ## backward sampling
        mm[t] = m[t] + C[t] / R[t+1] * (mm[t+1] - a[t+1])
        CC[t] = C[t] - (C[t]^2) / (R[t+1]^2) * (R[t+1] - CC[t+1])
        x[t]  = rnorm(1, m[t] + B[t] * (x[t+1] - a[t+1]), sqrt(H[t]))
    }
    return(list(x = x, m = m, C = C, mm = mm, CC = CC, llike = llike))
}

## Tobit model, parameter estimation by MCMC using FFBS
tobit_ffbs_mcmc <- function(y, c_1, ## data
                            a, b, c, d, m0, C0, ## priors
                            burn, M, step, ## arguments for mcmc
                            V1, W1, Phi1 ## initial value
                            ){
    ## functions, storages,
    niter = burn + M
    draws = NULL
    all_draws = NULL
    cs.ind <- which(y == c_1)
    ## iterations
    for (iter in 1:niter){
        run <- ffbs(y, V1, W1, Phi1, m0, C0)
        x <- run$x
        V1 <- 1 / rgamma(1, a + n / 2, b + sum((y - x)^2) / 2)
        W1 <- 1 / rgamma(1, c + (n - 1) / 2, d + sum((x[-1] - Phi1 * x[-length(x)])^2) / 2)
        regx <- lm(x[-length(x)] ~ 0 + x[-1])
        Phi1 = rnorm(1, coef(summary(regx))[1], coef(summary(regx))[2])
        ## impute censored value
        ## y[cs.ind] <- x[cs.ind] + rnorm(length(cs.ind), mean = 0, sd = sqrt(V1)) ## truncated normal
        ## y[cs.ind] <-  rtmvnorm(
        ##     n = 1,
        ##     mean = x[cs.ind],
        ##     sigma = sqrt(V1) * diag(length(cs.ind)),
        ##     upper = rep(c_1, length(cs.ind))
        ## )
        if (length(cs.ind) > 0){
            SampleY <- function(x) rtruncnorm(1, a = -Inf, b = c_1, mean = x, sd = V1)
            y[cs.ind] <- sapply(x[cs.ind], SampleY)
        }
        if (iter %% step == 0) {
            draws <- rbind(
                draws,
                c(V1, W1, Phi1, x, y)
            )
            cat(iter, '\n', 'V1: ', V1, "W1", W1, "Phi1", Phi1,
                ## "imputed y;", y[cs.ind],
                '\n')
        }
    }
    ## summary
    all_draws = draws[,1:3]
    draws = draws[(burn+1):(niter),]
    xs = draws[,4:(n+3)]
    ys = draws[, (n+4):(2 * n + 3)]
    Vs = draws[, 1]
    Ws = draws[, 2]
    Phis = draws[, 3]
    return(list(
        "xs" = xs,
        "ys" = ys,
        "Vs" = Vs,
        "Ws" = Ws,
        "Phis" = Phis
    ))
}

# Simulate states and data
W = 0.5
V = 1.0
n = 1000
m0 = 0.0 # mean of x_0
C0 = 10.0 # variance of x_0
x0 = 0
Phi = 0.8
## left bound for tobit model
## c_1 = -1 ## variance term gets larger and larger
c_1 = -2

trial <- function(seed){
    ## seed
    set.seed(seed)
    ## generate data
    w  = rnorm(n,0,sqrt(W))
    v  = rnorm(n,0,sqrt(V))
    x  = ystar = y = rep(0,n)
    x[1] = x0 + w[1]
    ystar[1] = x[1] + v[1] # observation before censoring.
    y[1] = max(ystar[1], c_1)
    for (t in 2:n){
        x[t] = Phi * x[t-1] + w[t]
        ystar[t] = x[t] + v[t]
        ## make it tobit
        y[t] = max(ystar[t], c_1)
    }
    ## censoring rate
    cs_rate <- sum(y == c_1) / n * 100
    ## Hyperparameters
    a = 0.01; b = 0.01; c = 0.01; d = 0.01 # MCMC step
    burn <- 100
    M <- 1000
    step <- 1
    V1 <- V ## initial value is true value
    W1 <- W ## initial value is true value
    Phi1 <- Phi
    ## mcmc_tobit
    sol <- tobit_ffbs_mcmc(y, c_1, ## data
                           a, b, c, d, m0, C0, ## priors
                           burn, M, step, ## arguments for mcmc
                           V1, W1, Phi1 ## initial value
                           )
    ## lx <- apply(sol$xs, 2, q025)ppp
    ## mx <- apply(sol$xs, 2, mean)
    ## ux <- apply(sol$xs, 2, q975)
    ## ly <- apply(sol$ys, 2, q025)
    ## my <- apply(sol$ys, 2, mean)
    ## uy <- apply(sol$ys, 2, q975)
    V <- mean(sol$Vs, na.rm = TRUE)
    W <- mean(sol$Ws, na.rm = TRUE)
    Phi <- mean(sol$Phis, na.rm = TRUE)
    ## mcmc_complete
    sol_comp <- tobit_ffbs_mcmc(ystar, -Inf, a, b, c, d, m0, C0, burn, M, step, V1, W1, Phi1)
    V_comp <- mean(sol_comp$Vs, na.rm = TRUE)
    W_comp <- mean(sol_comp$Ws, na.rm = TRUE)
    Phi_comp <- mean(sol_comp$Phis, na.rm = TRUE)
    return(c(cs_rate, V, W, Phi, V_comp, W_comp, Phi_comp))
}

## simulation
set.seed(123457)
nsamp <- 20
fwrite(data.table(t(c("seed", "cs_rate", "V", "W", "Phi", "V_comp", "W_comp", "Phi__comp"))), "simu_result.csv", append = TRUE)
for (seed in .Random.seed[1:nsamp]){
    cat(seed, '\n')
    simu <- trial(seed)
    rs <- c(seed, simu)
    fwrite(data.table(t(rs)), "simu_result.csv", append = TRUE)
}

## set.seed(20190601);
## set.seed(6022019)
## ## set.seed(123457)
## w  = rnorm(n,0,sqrt(W))
## v  = rnorm(n,0,sqrt(V))
## x  = ystar = y = rep(0,n)
## x[1] = x0 + w[1]
## ystar[1] = x[1] + v[1] # observation before censoring.
## y[1] = max(ystar[1], c_1)
## for (t in 2:n){
##     x[t] = Phi * x[t-1] + w[t]
##     ystar[t] = x[t] + v[t]
##     ## make it tobit
##     y[t] = max(ystar[t], c_1)
## }
## plot(x)
## plot(y)
## cat('left censoring rate: ', sum(y == c_1) / n * 100, '%\n')

## # actual smoother (for plotting)
## ks = Ksmooth0(num=n, ystar, A=1, m0, C0, Phi = Phi, cQ = sqrt(W), cR = sqrt(V))
## xsmooth = as.vector(ks$xs)

## ## true empirical contour graph
## run = ffbs(ystar, V, W, Phi, m0, C0)
## m = run$m;
## C = run$C;
## mm = run$mm
## CC = run$CC;
## L1 = m - 2 * C;
## U1 = m + 2 * C
## L2 = mm - 2 * CC;
## U2 = mm + 2 * CC
## N = 50
## Vs = seq(0.1, 2, length = N)
## Ws = seq(0.1, 2, length = N)
## likes = matrix(0,N,N)
## for (i in 1:N){
##     for (j in 1:N){
##         V = Vs[i]
##         W = Ws[j]
##         run = ffbs(ystar, V, W, Phi, m0, C0)
##         likes[i,j] = run$llike
##     }
## }

## ## Hyperparameters
## a = 0.01; b = 0.01; c = 0.01; d = 0.01 # MCMC step
## set.seed(90210)
## burn <- 10
## M <- 1000
## step <- 1
## V1 <- V ## initial value is true value
## W1 <- W ## initial value is true value
## Phi1 <- Phi



## test <- tobit_ffbs_mcmc(y, c_1, a, b, c, d, m0, C0, burn, M, step, V1, W1, Phi1)

## q025 = function(x){quantile(x,0.025)}
## q975 = function(x){quantile(x,0.975)}
## lx = apply(xs, 2, q025)
## mx = apply(xs, 2, mean)
## ux = apply(xs, 2, q975)
## ly = apply(ys, 2, q025)
## my = apply(ys, 2, mean)
## uy = apply(ys, 2, q975)

## ## plot of the data
## par(mfrow = c(2,2), mgp = c(1.6,.6,0), mar = c(3,3.2,1,1))
## ts.plot(ts(x), ts(y), ylab='', col = c(1, 8), lwd = 2)
## ts.plot(ts(x), ts(y), ts(ystar), ylab='', col = c(1, 8, 3), lwd = 2)
## points(ystar)
## legend(0, 11, legend = c("x(t)", "y(t)"), lty = 1, col = c(1, 8), lwd = 2, bty = "n",
##        pch = c(-1,1))
## contour(Vs, Ws, exp(likes), xlab = expression(sigma[v]^2),
##         ylab = expression(sigma[w]^2), drawlabels = FALSE, ylim = c(0,1.2))
## points(draws[,1:2], pch=16, col = rgb(.9,0,0,0.3), cex = .7)

## hist(draws[, 1], ylab = "Density", main = "", xlab = expression(sigma[v]^2))
## abline(v = mean(draws[, 1]), col=3, lwd=3)
## hist(draws[, 2], main = "", ylab = "Density", xlab = expression(sigma[w]^2))
## abline(v = mean(draws[, 2]), col=3, lwd=3)
## hist(draws[, 3], ylab = "Density", main = "", xlab = expression(Phi))
## abline(v = mean(draws[, 3], col = 3, lwd = 3))
## ## plot states

## par(mfrow = c(2,2))
## ts.plot(ts(draws[,1]))
## ts.plot(ts(draws[,2]))
## ts.plot(ts(draws[,3]))

## par(mgp=c(1.6,.6,0), mar=c(2,1,.5,0)+.5)
## plot(ts(mx), ylab='', type='n', ylim=c(min(y),max(y)))
## grid(lty=2);
## points(y)
## lines(xsmooth, lwd=4, col=rgb(1,0,1,alpha=.4))
## lines(mx, col= 4)
## xx=c(1:100, 100:1)
## yy=c(lx, rev(ux))
## polygon(xx, yy, border=NA, col= gray(.6,alpha=.2))
## lines(y, col=gray(.4))
## ## legend('topleft', c('true smoother', 'data', 'posterior mean', '95% of draws'),
## ##        lty=1, lwd=c(3,1,1,10), pch=c(-1,1,-1,-1), col=c(6, gray(.4), 4, gray(.6, alpha=.5)),
## ##        bg='white' )
