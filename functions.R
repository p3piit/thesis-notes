# This is the R script containing all the functions used for the thesis, the only purpose 
# of this script is to be sourced in all the Rmd files where the processes are explained 
# in a complite and detailed way. (Also has all the libraries)

# Libraries (probably I should use renv but I'm lazy)
library(tidyverse)
library(rpart)   
library(MASS)    
library(lme4)   
library(parallel) 
library(ggplot2) 

##############################################################
##############################################################

##############################################################
# From gmart_algo.Rmd
##############################################################

# b_fun: updates random effects b_i given current D and provided V_i
#        w is a LIST with one diagonal weight matrix W_i per cluster
b_fun <- function(G,    # number of clusters
                  Z,    # N x q random-effects design (stacked Z_i)
                  Vi,   # list length G with V_i matrices (n_i x n_i) 
                  D,    # q x q covariance of random effects
                  idx,  # list length G with row indices per cluster
                  y_t,  # length N pseudo-response (y_tilde)
                  fhat, # length N fixed-part prediction
                  w     # length N weights 
){
  b <- matrix(0, G, ncol(Z))
  for (g in seq_len(G)) {
    Zi   <- Z[idx[[g]], , drop = FALSE]
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                              # W_i^{1/2}
    
    Zi_w <- Wi_sq %*% Zi                                          # W^{1/2} Z_i
    rhs  <-  Wi_sq %*% y_t[idx[[g]]] - Wi_sq %*% fhat[idx[[g]]]   # W^{1/2} (y - fhat)
    
    b[g, ] <- D %*% t(Zi_w) %*% solve(Vi[[g]]) %*% rhs             # b_i = D (Z_w)^T V_i^{-1} rhs
  }
  b
}

##############################################################

# Vi_fun: builds list of V_i = Z_i D Z_i^T + sigma^2 I_{n_i}
Vi_fun <- function(Z,   # Z   : N x q random-effects design
                   D,   # D   : q x q covariance of random effects
                   s2,  # s2  : scalar sigma^2 (residual variance; homoskedastic)
                   G,   # G   : number of clusters
                   idx, # idx : list of row indices per cluster
                   w    # length N weights 
){
  Vi <- list()
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]                                    # Z_i
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                                     # W_i^{1/2}
    Zi_w <- Wi_sq %*% Zi                                                 # W^{1/2} Z_i
    Vi[[g]] <- Zi_w %*% D %*% t(Zi_w) + as.numeric(s2) * diag(nrow(Zi_w))  # V_i
  }
  # Output:
  # list of length G with n_i x n_i covariance matrices V_i
  return(Vi)
}

##############################################################

# sigma_fun: updates sigma^2 using paper's closed-form expression
sigma_fun <- function(N,    # N   : total number of observations
                      G,    # G   : number of clusters
                      idx,  # idx : list of row indices per cluster
                      b,    # b   : matrix q × n vector of random effects 
                      y,    # y   : response vector length N
                      Z,    # Z   : N x q random-effects design
                      D,    # D   : q x q covariance of random effects
                      Vi,   # Vi  : list of length G with V_i matrices (n_i x n_i)
                      s2,   # s2  : scalar sigma^2 (residual variance; homoskedastic)
                      fhat, # fhat: current fixed-effect fit
                      w     # length N weights 
){
  term_sum <- 0
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                    # W_i^{1/2}                         
    e <- Wi_sq %*% y[idx[[g]]] - Wi_sq %*% fhat[idx[[g]]] - Wi_sq %*% Zi%*% b[g, ]       #  weighted residuals epsilon_i = W_i^{1/2} y - W_i^{1/2} fhat - W_i^{1/2} Z_i b_i
    term_sum <- term_sum + crossprod(e) + s2 * (length(idx[[g]]) -
                                                  s2 * sum(diag(solve(Vi[[g]]))))     # add cluster-i contribution
  }
  # Output:
  #   scalar sigma^2 new estimate
  term_sum / N                                             # divide by total N
}

##############################################################

# D_fun: updates D using the paper's expression (average of b_i b_i^T plus shrinkage term)
D_fun <- function(G,   # G   : number of clusters
                  idx, # idx : list of row indices per cluster
                  b,   # b   : matrix q × n vector of random effects 
                  Z,   # Z   : N x q random-effects design
                  D,   # D   : q x q covariance of random effects
                  Vi,  # Vi  : list of length G with V_i matrices (n_i x n_i)
                  w    # length N weights 
){
  q <- ncol(Z)
  S <- matrix(0, q, q)                                     # accumulator
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]                      # Z_i
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                       # W_i^{1/2}
    Zi_w <- Wi_sq %*% Zi                                   # W^{1/2} Z_i
    S <- S + tcrossprod(b[g, ]) +                          # b_i b_i^T
      (D - D %*% t(Zi_w) %*% solve(Vi[[g]]) %*% Zi_w %*% D)    # shrinkage term per cluster
  }
  # Output:
  #   q x q updated D
  S/G                                                      # average over clusters
}

##############################################################

# gll_fun: computes generalized log-likelihood (up to additive constant)
gll_fun <- function(idx, # idx : list of row indices per cluster
                    b,   # b   : matrix q × n vector of random effects 
                    y,   # y   : response vector length N
                    Z,   # Z   : N x q random-effects design
                    D,   # D   : q x q covariance of random effects
                    s2,  # s2  : scalar sigma^2 (residual variance; homoskedastic)
                    fhat,# fhat: current fixed-effect fit
                    w    # length N weights 
){
  logdetD <- as.numeric(determinant(D, logarithm = TRUE)$modulus)     # log|D|
  term_b <- sum(rowSums((b %*% solve(D)) * b))                         # sum_i b_i^T D^{-1} b_i
  term_r <- 0
  term_logRi <- 0
  for (i in seq_along(idx)) {
    Zi <- Z[idx[[i]], , drop = FALSE]                      # Z_i
    Wi_sq <- diag(sqrt(w[idx[[i]]]))                       # W_i^{1/2}                         
    e <- Wi_sq %*% y[idx[[i]]] - Wi_sq %*% fhat[idx[[i]]] - Wi_sq %*% Zi%*% b[i, ]       #  weighted residuals epsilon_i = W_i^{1/2} y - W_i^{1/2} fhat - W_i^{1/2} Z_i b_i
    term_r <- term_r + sum(e^2) / s2                                  # (1/s2) * ||ri||^2  
    term_logRi <- term_logRi + length(idx[[i]]) * log(s2)              # n_i * log(s2)
  }
  
  # Output:
  #   scalar GLL value (up to constant)
  # GLL = sum_i [ (y - f - Zb)^T R_i^{-1} (y - f - Zb) + log|R_i| ] +
  # + sum_i b_i^T D^{-1} b_i + G*log|D|
  term_r + term_logRi + term_b + length(idx) * logdetD                 # assemble GLL
}

##############################################################

# ait_stop: Aitken acceleration based stopping rule
ait_stop <- function(gll,  #gll : numeric vector with GLL history (needs last three values)
                     tol,  #tol : tolerance for relative error to asymptote
                     it    #it  : current iteration index 
){
  converged <- FALSE
  l1 <- gll[it + 2] - gll[it + 1]                                      # L_t
  l2 <- gll[it + 1] - gll[it]                                          # L_{t-1}
  a <- l1 / l2                                                         # Aitken ratio
  if (is.finite(a) && a > 0 && a < 1) {                                # only meaningful if 0<a<1
    L_inf <- gll[it + 1] + l1 / (1 - a)                                # extrapolate asymptote
    err <- abs(L_inf - gll[it + 2]) / (abs(L_inf) + 1e-12)             # relative error to asymptote
    if (err < tol) converged <- TRUE                                   # declare convergence
  }
  converged                                                            # return flag
}

##############################################################

fit_gmert    <- function(df,               # df: data.frame with columns
                         #   id   : cluster identifier (factor or integer)
                         #   y    : numeric response
                         #   x1,x2,x3 : numeric covariates
                         #   leaf : true generating leaf (optional, for diagnostics)
                         max_iter_inn = 1000,  # maximum number of EM iterations (inner loop)
                         max_iter_out = 1000,  # maximum number of PQL iterations (outer loop)
                         tol = 1e-6,           # convergence tolerance for both loops (Aitken or relative diff)
                         cp = 0.0,             # rpart complexity parameter (pruning threshold)
                         minsplit = 50,        # minimum number of obs required to attempt a split
                         minbucket = 20,       # minimum number of obs in any terminal node
                         maxdepth = 5,         # maximum tree depth
                         xval = 10             # number of cross-validation folds in rpart
) {
  
  # --- Start timer ---
  time_start <- proc.time()                     # records user, system, elapsed time
  
  # --- Basic setup ---
  N <- nrow(df)                                 # total number of observations
  G <- length(unique(df$id))                    # number of clusters
  idx_by_cluster <- split(seq_len(N), df$id)    # list: row indices grouped by cluster
  
  y <- df$y                                     # response vector
  Z <- cbind(1, df$x1)                          # random-effects design: intercept + slope on x1
  q <- ncol(Z)                                  # number of random effects (q = 2)
  
  # --- Initialization (Step 0) ---
  M <- 0L                                       # outer-loop counter
  mu <- ifelse(y == 1, 0.75, 0.25)              # initial conditional means
  y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # initial pseudo-response (PQL linearization)
  w <- mu * (1 - mu)                            # initial working weights
  sigma2 <- 1                                   # initial residual variance
  D <- diag(2)                                  # initial random-effects covariance (identity)
  b <- matrix(0, G, q)                          # initialize cluster random effects
  gll <- c(0, 0)                                # GLL storage (2 slots for Aitken acceleration)
  eta_old <- rep(0, N)                          # previous eta for outer-loop convergence check
  converged_in <- c()                           # inner-loop convergence flags (per outer iteration)
  converged_out <- FALSE                        # outer-loop convergence flag
  
  # --- Outer loop (PQL updates) ---
  repeat {
    m = 0                                       # reset inner-loop counter
    
    # --- Inner loop (EM-like iterations) ---
    repeat {
      m <- m + 1L
      
      # (1.i) Partial E-step: compute adjusted pseudo-response y_tilde* = y_tilde - Z b
      zb <- numeric(N)                          # cluster-specific random contributions
      for (g in seq_len(G)) {
        idx <- idx_by_cluster[[g]]              # indices for cluster g
        zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]  # Z_i b_i
      }
      y_star <- y_t - zb                        # adjusted response for tree fit
      
      # (1.ii) M-step: fit regression tree for fixed effects f(X)
      ctrl <- rpart.control(cp = cp, minsplit = minsplit, xval = xval,
                            minbucket = minbucket, maxdepth = maxdepth)
      Xdf <- data.frame(x1 = df$x1, x2 = df$x2, x3 = df$x3)
      tree <- rpart(y_star ~ x1 + x2 + x3,
                    data = cbind(y_star = y_star, Xdf),
                    weights = w, method = "anova", control = ctrl)
      fhat <- as.numeric(predict(tree, newdata = Xdf))
      
      # (1.iii) Update random effects b_i
      Vi <- Vi_fun(Z = Z, D = D, s2 = sigma2, G = G, idx = idx_by_cluster, w = w)  # build V_i per cluster
      b <- b_fun(G = G, Z = Z, Vi = Vi, D = D, idx = idx_by_cluster,
                 y_t = y_t, fhat = fhat, w = w)                                   # update b_i estimates
      
      # (2.i) Update sigma^2 (residual variance)
      sigma2 <- sigma_fun(N = N, G = G, idx = idx_by_cluster,
                          b = b, y = y_t, Z = Z, D = D,
                          Vi = Vi, s2 = sigma2, fhat = fhat, w = w)
      
      # (2.ii) Update D (random-effects covariance)
      D <- D_fun(G = G, idx = idx_by_cluster,
                 Z = Z, D = D, Vi = Vi, b = b, w = w)
      
      # --- Inner-loop convergence check (GLL stabilization) ---
      gll[m + 2] <- gll_fun(D = D, b = b, idx = idx_by_cluster,
                            Z = Z, y = y_t, fhat = fhat, s2 = sigma2, w = w)
      if (m > 1L) {
        rel <- abs(gll[m + 2] - gll[m + 1]) / (abs(gll[m + 1]) + 1e-12)
        if (rel < tol) { n_iter <- m; converged_in_t <- TRUE; break }
      }
      
      if (m >= max_iter_inn) {                  # max-iteration guard
        n_iter <- m
        converged_in_t <- FALSE
        # message(sprintf("WARNING: the EM algorithm did not converge in %d iterations.", max_iter_inn))
        break
      }
    }
    
    # --- Outer-loop update (PQL step) ---
    converged_in <- c(converged_in, converged_in_t)
    zb <- numeric(N)
    for (g in seq_len(G)) {
      idx <- idx_by_cluster[[g]]
      zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]
    }
    eta <- fhat + zb                            # recompute linear predictor
    
    # (Outer stopping rule – paper style)
    d_eta <- sqrt(mean((eta - eta_old)^2))      # RMS change of eta
    if (d_eta < tol) {
      converged_out <- TRUE
      break
    }
    
    # Update working quantities for next outer iteration
    eta_old <- eta
    mu <- exp(eta) / (1 + exp(eta))             # updated conditional means
    y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # new pseudo-response
    w <- mu * (1 - mu)                          # new working weights
    
    M <- M + 1L
    if (M >= max_iter_out) {                    # guard against outer non-convergence
      converged_out <- FALSE
      message(sprintf("WARNING: the PQL algorithm did not converge in %d iterations.", max_iter_out))
      break
    }
  }
  
  # --- Stop timer and compute elapsed time ---
  time_end <- proc.time()
  elapsed <- as.numeric((time_end - time_start)["elapsed"])   # total runtime in seconds
  
  # --- Return fitted components ---
  list(
    tree = tree,                   # fitted rpart tree (fixed-effects function f(X))
    b = b,                         # estimated random effects (G x q)
    D = D,                         # estimated random-effects covariance
    sigma2 = sigma2,               # estimated residual variance
    mu = mu,                       # conditional means
    converged_in = converged_in,   # convergence flags for inner loops
    converged_out = converged_out, # convergence flag for outer loop
    n_iter = n_iter,               # iterations performed (inner loop)
    train_ids = df$id,             # cluster identifiers in training set
    gll = gll,                     # GLL trace (for diagnostics)
    tol = tol,                      # convergence tolerance used
    time = elapsed                 # total runtime (seconds)
  )
}

##############################################################

predict_gmert <- function(fit, new_df, thr = 0.5  
) {
  # fixed part: regression tree predictions using predictors only
  fhat <- as.numeric(predict(fit$tree, new_df))        # tree fitted on y*; use on new data
  # construct random effects design for new data: [1, x1]
  Znew <- cbind(1, new_df$x1)
  # container for random effects contributions
  add <- numeric(nrow(new_df))
  # clusters seen during training (so we have estimated b_i)
  clus_fit <- sort(unique(fit$train_ids))
  # map each new cluster id to its index in clus_fit (NA if unseen)
  map <- match(new_df$id, clus_fit)
  seen <- !is.na(map)                                  # TRUE for rows belonging to seen clusters
  if (any(seen)) {
    # add random effect contribution: row-wise sum of [1, x1] * b_i
    add[seen] <- rowSums(Znew[seen, , drop = FALSE] *
                           fit$b[map[seen], ])
  }
  # total prediction = fixed part + random effects contribution
  eta   <- fhat + add
  p     <- 1 / (1 + exp(-eta))
  
  # return 0/1; 
  yhat  <- ifelse(p >= thr, 1, 0)
  yhat
}

##############################################################
##############################################################

##############################################################
# From Simulation_notes.Rmd


##############################################################

sim_data_gmert <- function(G = 50,          # number of clusters
                           n_i = 40,        # observations per cluster
                           beta0 = 0.5,     # fixed intercept
                           beta1 = 1.2,     # fixed effect for x1
                           beta2 = -0.8,    # fixed effect for x2
                           beta3 = 0.6,     # fixed effect for x3
                           sigma_b0 = 0.8,  # SD of random intercept
                           sigma_b1 = 0.5,  # SD of random slope for x1
                           rho = 0.2,       # correlation between b0 and b1
                           seed = 123) {    # random seed for reproducibility
  
  set.seed(seed)
  
  # covariance matrix of random effects (b0_i, b1_i)
  D <- matrix(c(sigma_b0^2, rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)
  
  # generate random effects for G clusters
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)
  
  # cluster identifiers
  id <- rep(1:G, each = n_i)
  
  # generate covariates: x1 (uniform), x2 (normal), x3 (binary)
  x1 <- runif(G * n_i, -2, 2)
  x2 <- rnorm(G * n_i, 0, 1)
  x3 <- rbinom(G * n_i, 1, 0.5)
  
  # linear predictor: fixed part + random part (b0_i + b1_i * x1)
  eta <- numeric(G * n_i)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * x2[idx] + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }
  
  # convert to probabilities via logistic link
  p <- 1 / (1 + exp(-eta))
  
  # binary outcome drawn from Bernoulli(p)
  y <- rbinom(G * n_i, 1, p)
  
  # return as data.frame ready for model fitting
  data.frame(
    id = factor(id),
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}

##############################################################

# ---------------------------------------------------------------
# Train/test split per cluster using the sim_data_gmert() output
# ---------------------------------------------------------------
split_gmert_data <- function(df, train_prop = 0.7, seed = 123) {
  set.seed(seed)
  
  train_idx <- integer(0)
  test_idx  <- integer(0)
  ids <- unique(df$id)
  
  for (g in ids) {
    idx_g <- which(df$id == g)
    n_train_g <- floor(length(idx_g) * train_prop)
    train_g <- sample(idx_g, n_train_g)
    test_g  <- setdiff(idx_g, train_g)
    train_idx <- c(train_idx, train_g)
    test_idx  <- c(test_idx, test_g)
  }
  
  train_df <- df[train_idx, ]
  test_df  <- df[test_idx, ]
  list(train = train_df, test = test_df)
}

##############################################################
##############################################################

##############################################################
# From optimization.Rmd

##############################################################
# Ainv_fun: precomputes A_i^{-1} (and optionally Z_i^T W_i Z_i) for all clusters
Ainv_fun <- function(G,        # G      : number of clusters
                     Z,        # Z      : N x q random-effects design
                     w,        # W      : list of weigths 
                     D,        # D      : q x q covariance of random effects (current)
                     sigma2,   # sigma2 : residual variance (current)
                     idx      # idx    : list length G with row indices for each cluster
){
  Ainv_list <- vector("list", G)
  
  Dinv <- solve(D)                                 # D^{-1} (q x q) once
  
  for (g in seq_len(G)) {
    Zi  <- Z[idx[[g]], , drop = FALSE]             # n_i x q
    Wi  <- diag(w[idx[[g]]])                      
    ZtWZ <- t(Zi) %*% Wi %*% Zi                     # Z_i^T W_i Z_i (q x q)
    
    A <- Dinv + (1 / sigma2) * ZtWZ                # A_i
    Ainv <- solve(A)                              # A_i^{-1} 
    
    Ainv_list[[g]] <- Ainv
  }
  
  Ainv_list
}

##############################################################

# b_fun_Ainv: updates random effects b_i using precomputed A_i^{-1}
b_fun_small <- function(G,        # G      : number of clusters
                        Z,        # Z      : N x q random-effects design
                        w,        # W      : list of weigths 
                        idx,      # idx    : list length G with row indices for each cluster
                        y_t,      # y_t    : pseudo-response vector (length N)
                        fhat,     # fhat   : fitted fixed-part prediction (length N)
                        Ainv,     # Ainv   : list length G; each Ainv[[g]] = A_i^{-1} (q x q)
                        sigma2    # sigma2 : residual variance (current)
){
  q <- ncol(Z)
  b <- matrix(0, G, q)
  
  for (g in seq_len(G)) {
    Zi   <- Z[idx[[g]], , drop = FALSE]           # n_i x q
    Wi   <- diag(w[idx[[g]]])                     # weights vector (length n_i)
    resid<- y_t[idx[[g]]] - fhat[idx[[g]]]        # (y_t - fhat) on cluster i
    
    rhs  <- (1 / sigma2) * t(Zi) %*% (Wi %*% resid)  # right-hand side
    b[g, ] <- as.vector(Ainv[[g]] %*% rhs)           # b_i = A_i^{-1} * rhs
  }
  
  b
}

##############################################################

# sigma_fun_Ainv: updates residual variance sigma2 using precomputed A_i^{-1}
sigma_fun_small <- function(N,        # N      : total number of rows
                            G,        # G      : number of clusters
                            idx,      # idx    : list length G with row indices for each cluster
                            Z,        # Z      : N x q random-effects design
                            w,        # W      : list of weigths
                            y_t,      # y_t    : pseudo-response vector
                            fhat,     # fhat   : fitted fixed-part prediction
                            b,        # b      : matrix G x q with current random effects
                            Ainv,     # Ainv   : list length G; each Ainv[[g]] = A_i^{-1}
                            sigma2   # sigma2 : residual variance (current)
){
  rss_total <- 0
  
  for (g in seq_len(G)) {
    Zi  <- Z[idx[[g]], , drop = FALSE]             # n_i x q
    wi  <- diag(w[idx[[g]]])                             # weights matrix 
    ni  <- length(wi)
    
    # residuals eps_i
    eps_i <- diag(sqrt(w[idx[[g]]])) %*% (y_t[idx[[g]]] - fhat[idx[[g]]] - as.vector(Zi %*% b[g, ]))
    
    # tr(V_i^{-1})
    trVi_inv <- (ni / sigma2) - (1 / sigma2^2) * sum(diag(Ainv[[g]] %*% t(Zi) %*% wi %*% Zi))
    
    # accumulate rss
    rss_total <- rss_total + t(eps_i) %% eps_i + sigma2 * (ni - sigma2 * trVi_inv)
  }
  
  rss_total / N
}

##############################################################

# D_fun_Ainv: updates random-effects covariance D using precomputed A_i^{-1}
D_fun_small <- function(G,        # G      : number of clusters
                        b,        # b      : matrix G x q with current random effects
                        Ainv      # Ainv   : list length G; each Ainv[[g]] = A_i^{-1}
){
  q <- ncol(b)
  D_new <- matrix(0, q, q)
  
  for (g in seq_len(G)) {
    D_new <- D_new + ( b %*% t(b) + Ainv[[g]] )   # b_i b_i^T + A_i^{-1}
  }
  
  D_new / G
}

##############################################################

# Main fitting function using Ainv-based updates
fit_gmert_small    <- function(df,         # df: data.frame with columns
                                                     # id   : cluster identifier (factor or integer)
                                                     # y    : numeric response
                                                     # x1,x2,x3 : numeric covariates
                                                     # leaf : true generating leaf (optional, for diagnostics)
                               max_iter_inn = 1000,  # maximum number of EM iterations (inner loop)
                               max_iter_out = 1000,  # maximum number of PQL iterations (outer loop)
                               tol = 1e-6,           # convergence tolerance for both loops (Aitken or relative diff)
                               cp = 0.0,             # rpart complexity parameter (pruning threshold)
                               minsplit = 50,        # minimum number of obs required to attempt a split
                               minbucket = 20,       # minimum number of obs in any terminal node
                               maxdepth = 5,         # maximum tree depth
                               xval = 10             # number of cross-validation folds in rpart
) {
  
  # --- Start timer ---
  time_start <- proc.time()              # records user, system, elapsed time
  
  # --- Basic setup ---
  N <- nrow(df)                                 # total number of observations
  G <- length(unique(df$id))                    # number of clusters
  idx_by_cluster <- split(seq_len(N), df$id)    # list: row indices grouped by cluster
  
  y <- df$y                                     # response vector
  Z <- cbind(1, df$x1)                          # random-effects design: intercept + slope on x1
  q <- ncol(Z)                                  # number of random effects (q = 2)
  
  # --- Initialization (Step 0) ---
  M <- 0L                                       # outer-loop counter
  mu <- ifelse(y == 1, 0.75, 0.25)              # initial conditional means
  y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # initial pseudo-response (PQL linearization)
  w <- mu * (1 - mu)                            # initial working weights
  sigma2 <- 1                                   # initial residual variance
  D <- diag(2)                                  # initial random-effects covariance (identity)
  b <- matrix(0, G, q)                          # initialize cluster random effects
  gll <- c(0, 0)                                # GLL storage (2 slots for Aitken acceleration)
  eta_old <- rep(0, N)                          # previous eta for outer-loop convergence check
  converged_in <- c()                           # inner-loop convergence flags (per outer iteration)
  converged_out <- FALSE                        # outer-loop convergence flag
  
  # --- Outer loop (PQL updates) ---
  repeat {
    m = 0                                       # reset inner-loop counter
    
    # --- Inner loop (EM-like iterations) ---
    repeat {
      m <- m + 1L
      
      # (1.i) Partial E-step: compute adjusted pseudo-response y_tilde* = y_tilde - Z b
      zb <- numeric(N)                          # cluster-specific random contributions
      for (g in seq_len(G)) {
        idx <- idx_by_cluster[[g]]              # indices for cluster g
        zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]  # Z_i b_i
      }
      y_star <- y_t - zb                        # adjusted response for tree fit
      
      # (1.ii) M-step: fit regression tree for fixed effects f(X)
      ctrl <- rpart.control(cp = cp, minsplit = minsplit, xval = xval,
                            minbucket = minbucket, maxdepth = maxdepth)
      Xdf <- data.frame(x1 = df$x1, x2 = df$x2, x3 = df$x3)
      tree <- rpart(y_star ~ x1 + x2 + x3,
                    data = cbind(y_star = y_star, Xdf),
                    weights = w, method = "anova", control = ctrl)
      fhat <- as.numeric(predict(tree, newdata = Xdf))
      
      # (1.iii) Update random effects b_i
      Ainv <- Ainv_fun(Z = Z, D = D, sigma2 = sigma2, G = G, idx = idx_by_cluster, w = w)  # build V_i per cluster
      b <- b_fun_small(G = G, Z = Z, idx = idx_by_cluster,
                       y_t = y_t, fhat = fhat, w = w, sigma2 = sigma2, Ainv = Ainv)                                   # update b_i estimates
      
      # (2.i) Update sigma^2 (residual variance)
      sigma2 <- sigma_fun_small(N = N, G = G, idx = idx_by_cluster,
                                b = b, y = y_t, Z = Z, 
                                Ainv = Ainv, sigma2 = sigma2, fhat = fhat, w = w)
      
      # (2.ii) Update D (random-effects covariance)
      D <- D_fun_small(G = G, b = b, Ainv = Ainv)
      
      # --- Inner-loop convergence check (GLL stabilization) ---
      gll[m + 2] <- gll_fun(D = D, b = b, idx = idx_by_cluster,
                            Z = Z, y = y_t, fhat = fhat, s2 = sigma2, w = w)
      if (m > 1L) {
        rel <- abs(gll[m + 2] - gll[m + 1]) / (abs(gll[m + 1]) + 1e-12)
        if (rel < tol) { n_iter <- m; converged_in_t <- TRUE; break }
      }
      
      if (m >= max_iter_inn) {                  # max-iteration guard
        n_iter <- m
        converged_in_t <- FALSE
        # message(sprintf("WARNING: the EM algorithm did not converge in %d iterations.", max_iter_inn))
        break
      }
    }
    
    # --- Outer-loop update (PQL step) ---
    converged_in <- c(converged_in, converged_in_t)
    zb <- numeric(N)
    for (g in seq_len(G)) {
      idx <- idx_by_cluster[[g]]
      zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]
    }
    eta <- fhat + zb                            # recompute linear predictor
    
    # (Outer stopping rule – paper style)
    d_eta <- sqrt(mean((eta - eta_old)^2))      # RMS change of eta
    if (d_eta < tol) {
      converged_out <- TRUE
      break
    }
    
    # Update working quantities for next outer iteration
    eta_old <- eta
    mu <- exp(eta) / (1 + exp(eta))             # updated conditional means
    y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # new pseudo-response
    w <- mu * (1 - mu)                          # new working weights
    
    M <- M + 1L
    if (M >= max_iter_out) {                    # guard against outer non-convergence
      converged_out <- FALSE
      message(sprintf("WARNING: the PQL algorithm did not converge in %d iterations.", max_iter_out))
      break
    }
  }
  
  # --- Stop timer and compute elapsed time ---
  time_end <- proc.time()
  elapsed <- as.numeric((time_end - time_start)["elapsed"])   # total runtime in seconds
  
  
  # --- Return fitted components ---
  list(
    tree = tree,                   # fitted rpart tree (fixed-effects function f(X))
    b = b,                         # estimated random effects (G x q)
    D = D,                         # estimated random-effects covariance
    sigma2 = sigma2,               # estimated residual variance
    mu = mu,                       # conditional means
    converged_in = converged_in,   # convergence flags for inner loops
    converged_out = converged_out, # convergence flag for outer loop
    n_iter = n_iter,               # iterations performed (inner loop)
    train_ids = df$id,             # cluster identifiers in training set
    gll = gll,                     # GLL trace (for diagnostics)
    tol = tol,                     # convergence tolerance used
    time = elapsed                 # total runtime (seconds)
  )
}







































