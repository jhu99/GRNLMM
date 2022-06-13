eigen2 <- function (spd, decreasing = FALSE) 
{
  foo <- eigen(spd)
  bar <- foo
  bar$values <- foo$values[order(foo$values, decreasing = decreasing)]
  bar$vectors <- foo$vectors[, order(foo$values, decreasing = decreasing)]
  return(bar)
}

eigen_p <- function (V_g, V_e, tol = 1/10000) 
{
  V_g<-Re(V_g)
  V_e<-Re(V_e)
  stopifnot(isSymmetric(V_g, tol = tol), isSymmetric(V_e, tol = tol))
  d_size <- nrow(V_g)
  logdet_Ve <- 0
  Lambda <- matrix(nrow = d_size, ncol = d_size)
  V_e_temp <- matrix(nrow = d_size, ncol = d_size)
  V_e_h <- matrix(0, nrow = d_size, ncol = d_size)
  V_e_hi <- matrix(0, nrow = d_size, ncol = d_size)
  VgVehi <- matrix(nrow = d_size, ncol = d_size)
  U_l <- matrix(nrow = d_size, ncol = d_size)
  V_e_temp <- V_e
  eout <- eigen2(V_e_temp)
  D_l <- eout$values
  U_l <- eout$vectors
  if (length(U_l == 1)) 
    U_l <- as.matrix(U_l)
  for (i in 1:d_size) {
    d <- D_l[i]
    if (d > 0) {
      logdet_Ve <- logdet_Ve + log(d)
      U_col <- U_l[, i]
      d <- sqrt(d)
      V_e_h <- V_e_h + d * U_col %*% t(U_col)
      V_e_hi <- V_e_hi + U_col %*% t(U_col)/d
    }
  }
  VgVehi <- V_g %*% V_e_hi
  Lambda <- V_e_hi %*% VgVehi
  eout <- eigen2(Lambda)
  D_l <- eout$values
  U_l <- eout$vectors
  if (length(U_l) == 1) 
    U_l <- as.matrix(U_l)
  D_l<-Re(D_l)
  D_l[D_l < 0] <- 0
  UltVeh <- t(U_l) %*% V_e_h
  UltVehi <- t(U_l) %*% V_e_hi
  return(list(logdet_Ve = logdet_Ve, UltVeh = UltVeh, UltVehi = UltVehi, 
              D_l = D_l))
}

calc_XHiY <- function (eval, D_l, X, UltVehiY) 
{
  stopifnot(length(eval) == ncol(X), is.vector(eval), is.vector(D_l))
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  xHiy <- rep(0, d_size * c_size)
  for (i in 1:d_size) {
    dl <- D_l[i]
    for (j in 1:c_size) {
      d <- 0
      for (k in 1:n_size) {
        x <- X[j, k]
        y <- UltVehiY[i, k]
        delta <- eval[k]
        d <- d + x * y/(delta * dl + 1)
      }
      xHiy[(j - 1) * d_size + i] <- d
    }
  }
  return(xHiy)
}

MphCalcLogL <- function (eval, D_l, Qi, UltVehiY, xHiy) 
{
  n_size <- length(eval)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  logl <- 0
  for (k in 1:n_size) {
    delta <- eval[k]
    for (i in 1:d_size) {
      y <- UltVehiY[i, k]
      dl <- D_l[i]
      d <- delta * dl + 1
      logl <- logl + y^2/d + log(d)
    }
  }
  Qiv <- Qi %*% xHiy
  d <- t(xHiy) %*% Qiv
  stopifnot(length(d) == 1)
  logl <- logl - d
  return(-0.5 * logl)
}

calc_omega <- function (eval, D_l) 
{
  n_size <- length(eval)
  d_size <- length(D_l)
  OmegaU <- matrix(nrow = d_size, ncol = n_size)
  OmegaE <- OmegaU
  for (k in 1:n_size) {
    delta <- eval[k]
    for (i in 1:d_size) {
      dl <- D_l[i]
      d_u <- dl/(delta * dl + 1)
      d_e <- d_u * delta
      OmegaU[i, k] <- d_u
      OmegaE[i, k] <- d_e
    }
  }
  return(list(OmegaU, OmegaE))
}

UpdateRL_B <- function (xHiy, Qi, d_size) 
{
  dc_size <- nrow(Qi)
  c_size <- dc_size/d_size
  b <- Qi %*% xHiy
  UltVehiB <- matrix(nrow = d_size, ncol = c_size)
  for (i in 1:c_size) {
    b_subcol <- b[(1 + (i - 1) * d_size):(i * d_size)]
    UltVehiB[, i] <- b_subcol
  }
  return(UltVehiB)
}

update_u <- function (OmegaE, UltVehiY, UltVehiBX) 
{
  UltVehiU <- UltVehiY
  UltVehiU <- UltVehiU - UltVehiBX
  UltVehiU <- UltVehiU * OmegaE
  return(UltVehiU)
}

update_e <- function (UltVehiY, UltVehiBX, UltVehiU) 
{
  UltVehiE <- UltVehiY - UltVehiBX - UltVehiU
  return(UltVehiE)
}

calc_sigma <- function (eval, D_l, X, OmegaU, OmegaE, UltVeh, Qi) 
{
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- length(D_l)
  dc_size <- nrow(Qi)
  Sigma_ee <- matrix(0, nrow = d_size, ncol = d_size)
  Sigma_uu <- Sigma_ee
  for (k in 1:n_size) {
    OmegaU_col <- OmegaU[, k]
    OmegaE_col <- OmegaE[, k]
    diag(Sigma_uu) <- diag(Sigma_uu) + OmegaU_col
    diag(Sigma_ee) <- diag(Sigma_ee) + OmegaE_col
  }
  M_u <- matrix(0, nrow = dc_size, ncol = d_size)
  M_e <- M_u
  for (k in 1:n_size) {
    delta <- eval[k]
    for (i in 1:d_size) {
      dl <- D_l[i]
      for (j in 1:c_size) {
        x <- X[j, k]
        d <- x/(delta * dl + 1)
        M_e[(j - 1) * d_size + i, i] <- d
        M_u[(j - 1) * d_size + i, i] <- d * dl
      }
    }
    QiM <- Qi %*% M_u
    Sigma_uu <- Sigma_uu + t(M_u) %*% QiM * delta
    QiM <- Qi %*% M_e
    Sigma_ee <- Sigma_ee + t(M_e) %*% QiM
  }
  M <- Sigma_uu %*% UltVeh
  Sigma_uu <- t(UltVeh) %*% M
  M <- Sigma_ee %*% UltVeh
  Sigma_ee <- t(UltVeh) %*% M
  return(list(Sigma_ee, Sigma_uu))
}

update_v <- function (eval, U, E, Sigma_uu, Sigma_ee, tol = 1/10000) 
{
  stopifnot(isSymmetric(Sigma_uu, tol = tol), isSymmetric(Sigma_ee, 
                                                          tol = tol))
  n_size <- length(eval)
  d_size <- nrow(U)
  V_g <- matrix(0, nrow = d_size, ncol = d_size)
  V_e <- V_g
  for (k in 1:n_size) {
    delta <- eval[k]
    U_col <- U[, k]
    V_g <- V_g + U_col %*% t(U_col)/delta
  }
  V_e <- E %*% t(E)
  V_g <- V_g + Sigma_uu
  V_e <- V_e + Sigma_ee
  V_g <- V_g/n_size
  V_e <- V_e/n_size
  return(list(V_g, V_e))
}

grnlmm <- function(x,V_g,V_e)
{
  K <- cor(x)
  eig <- eigen2(K)
  M<-rowMeans(x)
  m<-matrix(nrow = nrow(x),ncol = ncol(x))
  for (i in 1:nrow(m)) {
    m[i,]<-M[i]
  }
  Y<-x-m
  X<-matrix(0,1,ncol = ncol(x))
  rm(m)
  rm(x)
  Y<-as.matrix(Y)
  max_iter <- 10000
  max_prec <- 1/1e+06
  eval <- eig$values
  verbose_output <- FALSE
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- nrow(Y)
  dc_size <- c_size * d_size
  XXt<-matrix(1,1,1)
  logl_const <- -(n_size - c_size) * d_size * log(2 * pi)/2 + 
    d_size * log(det(XXt))/2
  out <- list()
  for (t in 1:max_iter) {
    ep_out <- eigen_p(V_g, V_e)
    logdet_Ve <- ep_out[[1]]
    UltVeh <- ep_out[[2]]
    UltVehi <- ep_out[[3]]
    D_l <- ep_out[[4]]
    Qi <-matrix(0,nrow(V_g),nrow(V_g))
    lndetQ <- 1
    UltVehiY <- UltVehi %*% Y
    xHiy <- calc_XHiY(eval, D_l, X, UltVehiY)
    logl_new <- logl_const + MphCalcLogL(eval = eval, xHiy = xHiy, 
                                                  D_l = D_l, UltVehiY = UltVehiY, Qi = Qi) - 0.5 * 
      n_size * logdet_Ve
    logl_new <- logl_new - 0.5 * (lndetQ - c_size * logdet_Ve)
    if (t > 1) {
      logl_new<-Re(logl_new)
      if (logl_new - logl_old < max_prec) {
        break
      }
    }
    logl_old <- logl_new
    co_out <- calc_omega(eval, D_l)
    OmegaU <- co_out[[1]]
    OmegaE <- co_out[[2]]
    UltVehiB <- UpdateRL_B(xHiy, Qi, d_size = nrow(Y))
    UltVehiBX <- UltVehiB %*% X
    UltVehiU <- update_u(OmegaE, UltVehiY, UltVehiBX)
    UltVehiE <- update_e(UltVehiY, UltVehiBX, UltVehiU)
    U_hat <- t(UltVeh) %*% UltVehiU
    E_hat <- t(UltVeh) %*% UltVehiE
    B <- t(UltVeh) %*% UltVehiB
    cs_out <- calc_sigma(eval = eval, D_l = D_l, X = X, OmegaU = OmegaU, 
                                  OmegaE = OmegaE, UltVeh = UltVeh, Qi = Qi)
    Sigma_ee <- cs_out[[1]]
    Sigma_uu <- cs_out[[2]]
    Sigma_ee <- Re(Sigma_ee)
    Sigma_uu <- Re(Sigma_uu)
    uv_out <- update_v(eval, U_hat, E_hat, Sigma_uu, Sigma_ee)
    V_g <- uv_out[[1]]
    V_e <- uv_out[[2]]
    if (verbose_output) {
      out[[t]] <- list(logl_new = logl_new, Vg = V_g, Ve = V_e, 
                       Sigma_uu = Sigma_uu, Sigma_ee = Sigma_ee, B = B, 
                       U_hat = U_hat, E_hat = E_hat, OmegaU = OmegaU, 
                       OmegaE = OmegaE, logdet_Ve = logdet_Ve, UltVeh = UltVeh, 
                       UltVehi = UltVehi, Dl = D_l, xHiy = xHiy, logl_const = logl_const, 
                       UltVehiU = UltVehiU)
    }
    else {
      out[[1]] <- list(logl_new = logl_new, Vg = V_g, Ve = V_e, 
                       Sigma_uu = Sigma_uu, Sigma_ee = Sigma_ee, B = B, 
                       U_hat = U_hat, E_hat = E_hat, OmegaU = OmegaU, 
                       OmegaE = OmegaE, logdet_Ve = logdet_Ve, UltVeh = UltVeh, 
                       UltVehi = UltVehi, Dl = D_l, xHiy = xHiy, logl_const = logl_const, 
                       UltVehiU = UltVehiU)
    }
  }
  if (length(out) == max_iter) {
    warning("EM algorithm didn't converge.")
  }
  return(out[[1]][["Vg"]])
}
x<-as.matrix(read.table("data.txt", sep="\t"))
g<-as.matrix(read.table("tf.txt",sep = "\t"))
row.names(x)<-g[,1]
ma<-matrix(0,nrow = nrow(x),ncol = nrow(x))
for (i in 1:nrow(x)) {
  ma[i,i]<-1
}
ma<-as.matrix(ma)
V_g <- ma
V_e <- ma
zz<-grnlmm(x,V_g,V_e)
