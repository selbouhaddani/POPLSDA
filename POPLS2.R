POPLS.Estep2 <- function(data, prm){

  if(!is.list(data)) stop("need a correct data object")
  if(mean(c("W","P","beta","SigT","SigU","sig2e","sig2eps") %in% names(prm)) != 1) stop("need a correct prm object")

  p <- nrow(prm$W)
  r <- ncol(prm$W)
  rx <- ncol(prm$P[[1]]) #* sign(ssq(prm$P[[1]])) #function(Ps) {ncol(Ps)*sign(ssq(Ps))}
  #prm$P %<>% lapply(`*`, ssq(prm$P))

  ### old
  #sig2e_inv <- diag(1/unlist(prm$sig2e), length(prm$sig2e))
  #sig2eps_inv <- diag(1/unlist(prm$sig2eps), length(prm$sig2eps))
  var_T <- sapply(simplify = FALSE, 1:length(data),
                  function(ii) {
                    solve(with(prm, diag(1,r) +
                                 1/sig2e[[ii]]*crossprod(W) +
                                 1/sig2eps[[ii]]*tcrossprod(beta)))
                  }) #%>% Reduce(f=`+`) %>% divide_by(length(data))

  mu_T <- sapply(simplify = FALSE, 1:length(data),
                 function(ii) {
                   data[[ii]]$X%*%(prm$W/prm$sig2e[[ii]])%*%var_T[[ii]] +
                     as.matrix(data[[ii]]$Y) %*%
                     t(prm$beta/prm$sig2eps[[ii]])%*%var_T[[ii]]
                 })

  Stt <- sapply(simplify = FALSE, 1:length(mu_T), function(ii) nrow(mu_T[[ii]])*var_T[[ii]] + crossprod(mu_T[[ii]]))
  ###

  var_Z <- sapply(simplify = FALSE, 1:length(data), function(ii) {
    Ps <- prm$P[[ii]]
    sig2e_s <- prm$sig2e[[ii]]
    sig2eps_s <- prm$sig2eps[[ii]]
    Gamma <- with(prm, cbind(rbind(W, t(beta)),rbind(Ps, matrix(0,1,ncol(Ps)))) )
    GammaTe <- with(prm, cbind(rbind(1/sig2e_s*W, 1/sig2eps_s*t(beta)),rbind(1/sig2e_s*Ps, matrix(0,1,ncol(Ps)))) )
    B = with(prm, crossprod(Gamma, GammaTe) )
    return(solve(diag(1, r+rx) + B))
  })

  mu_Z <- sapply(simplify = FALSE, 1:length(data), function(ii) {
    e <- data[[ii]]
    Ps <- prm$P[[ii]]
    sig2e_s <- prm$sig2e[[ii]]
    sig2eps_s <- prm$sig2eps[[ii]]
    Gamma <- with(prm, cbind(rbind(W, t(beta)),rbind(Ps, matrix(0,1,ncol(Ps)))) )
    cbind(e$X/sig2e_s, e$Y/sig2eps_s) %*% Gamma %*% var_Z[[ii]]
  })

  Szz <- sapply(simplify = FALSE, 1:length(data), function(ii)
    nrow(mu_Z[[ii]])*var_Z[[ii]] + crossprod(mu_Z[[ii]]))

  outp <- list(
    mu_T = sapply(simplify=F, mu_Z, function(e) e[,1:r] %>% as.matrix),
    mu_U = sapply(simplify=F, mu_Z, function(e) e[,-(1:r)] %>% as.matrix),
    var_T = sapply(simplify=F, var_Z, function(e) e[1:r,1:r] %>% as.matrix),
    var_U = sapply(simplify=F, var_Z, function(e) e[-(1:r),-(1:r)] %>% as.matrix),
    var_UT = sapply(simplify=F, var_Z, function(e) e[-(1:r),1:r] %>% as.matrix),
    Stt = sapply(simplify=F, Szz, function(e) e[1:r,1:r] %>% as.matrix),
    Suu = sapply(simplify=F, Szz, function(e) e[-(1:r),-(1:r)] %>% as.matrix),
    Sut = sapply(simplify=F, Szz, function(e) e[-(1:r),1:r] %>% as.matrix)
  )

  #outp$var_T <- Reduce(outp$var_T, f=`+`) %>% divide_by(length(data))
  var_T <- outp$var_T
  mu_T <- outp$mu_T

  See_trace <- sapply(simplify = F, 1:length(data), function(ii){
    with(data[[ii]], with(prm, {
      nrow(X)*tr(t(W)%*%W%*%var_T[[ii]]) +
        ( ssq(X) - 2*tr(t(X%*%W)%*%mu_T[[ii]]) + tr(t(W)%*%W%*%crossprod(mu_T[[ii]])) )
    }))
  })

  Sepseps <- sapply(simplify = F, 1:length(data), function(ii){
    with(data[[ii]], with(prm, {
      nrow(X)*t(beta)%*%var_T[[ii]]%*%beta +
        crossprod(Y - mu_T[[ii]]%*%beta)
    }))
  })

  loglik <- sapply(1:length(data), function(i){
    X <- data[[i]]$X
    Y <- data[[i]]$Y
    N <- nrow(X)
    p <- ncol(X)
    Ps <- prm$P[[i]]
    sig2e_s <- prm$sig2e[[i]]
    sig2eps_s <- prm$sig2eps[[i]]

    Gamma <- with(prm, cbind(rbind(W, t(beta)),rbind(Ps, matrix(0,1,ncol(Ps)))) );

    det.2 <- with(prm, c(determinant( diag(1,r+ncol(Ps)) +
                                        blockm(1/sig2e_s*crossprod(W)+1/sig2eps_s*tcrossprod(beta),
                                               1/sig2e_s*crossprod(W,Ps), 1/sig2e_s*crossprod(Ps)) )$mod) +
                    p*log(sig2e_s) + log(sig2eps_s) )
    tr.2 <- with(prm, ssq(cbind(X/sqrt(sig2e_s), Y/sqrt(sig2eps_s))) -
                   tr( crossprod(cbind(X/sig2e_s,Y/sig2eps_s)%*%Gamma) %*%var_Z[[i]] ) )

    logl.2 <- -0.5*(N*(p+1)*log(2*pi) + N*det.2 + tr.2)

    return(logl.2)
  })


  return(c(outp, list(See_trace=See_trace, Sepseps=Sepseps, logl = loglik)))
}



POPLS.Mstep <- function(Estep, data, prm){
  if(!is.list(data)) stop("need a correct data object")
  N_s <- sapply(data, function(e) nrow(e$X))
  sumN <- sum(N_s)
  outp <- prm
  rx <- ncol(outp$P[[1]]) #* sign(ssq(outp$P[[1]])) #function(Ps) {ncol(Ps)*sign(ssq(Ps))}

  XtT <- lapply(1:length(data),
                function(i) crossprod(data[[i]]$X, Estep$mu_T[[i]]) -
                  outp$P[[i]]%*%crossprod(Estep$mu_U[[i]],Estep$mu_T[[i]])) %>% Reduce(f="+")
  TtT <- lapply(1:length(data), function(i) Estep$Stt[[i]] ) %>% Reduce(f="+")
  What <- XtT %*% solve(TtT)

  XtU <- lapply(1:length(data),
                function(i) crossprod(data[[i]]$X, Estep$mu_U[[i]]) -
                  outp$W%*%crossprod(Estep$mu_T[[i]],Estep$mu_U[[i]]) )
  UtU <- lapply(1:length(data), function(i) Estep$Suu[[i]] )
  Phat <- lapply(1:length(data), function(i) XtU[[i]] %*% solve(UtU[[i]]))

  outp$W <- What#with(svd(What), u %*% diag(d,length(d)))
  outp$P <- lapply(Phat, function(e) e * sign(rx))

  YtT <- lapply(1:length(data),
                function(i) crossprod(data[[i]]$Y, Estep$mu_T[[i]]) ) %>% Reduce(f="+")
  outp$beta <- t(YtT %*% solve(TtT))

  sig2e_hat <- Estep$See_trace #Reduce(Estep$See_trace, f = `+`)
  outp$sig2e <- sapply(1:length(sig2e_hat),
                       function(ii) {sig2e_hat[[ii]] / nrow(What) / N_s[[ii]]})
  sig2eps_hat <- Estep$Sepseps #Reduce(Estep$Sepseps, f = `+`) %>% c %>% divide_by(sumN)
  outp$sig2eps <- sapply(1:length(sig2eps_hat),
                      function(ii) {sig2eps_hat[[ii]] / N_s[[ii]]})

  #sig2e_hat <- Reduce(Estep$See_trace, f = `+`)
  #outp$sig2e <- rep(sig2e_hat / nrow(What) / sumN , length(data))
  # outp$sig2eps <- Reduce(Estep$Sepseps, f = `+`) %>% c %>% divide_by(sumN) %>%
  #   rep(length(data))

  # outp$SigT <- Reduce(Estep$Stt, f = "+")
  # outp$SigT <- outp$SigT %>% divide_by(nrow(Estep$Stt[[1]])*sumN) %>%
  #   multiply_by(diag(1,nrow(Estep$Stt[[1]])))

  return(list(param = outp, mu_T=Estep$mu_T))
}



POPLS.EM <- function(data, prm, nsteps, atol){
  #prm.cur$SigT %<>% multiply_by(2)
  prm.cur <- prm
  logl <- matrix(NA, nsteps, length(data))
  logl2 <- NA*1:nsteps
  for(iter_i in 1:nsteps){
    #cat(iter_i,"--")
    E_i <- POPLS.Estep2(data, prm.cur)
    M_i <- POPLS.Mstep(E_i, data, prm.cur)
    prm.cur <- M_i$param
    prm.cur$W %<>% sweep(2,sign(prm.cur$beta),FUN = `*`)
    prm.cur$beta %<>% abs
    logl[iter_i,] <- (E_i$logl)
    logl2[iter_i] <- sum(E_i$logl)

    if(iter_i > 1 && abs(logl2[iter_i] - logl2[iter_i-1]) < atol) break
    #cat(iter_i, ", ", round(sum(E_i$logl)), "  |  ", sep="")
    #on.exit(return(list(param = prm.cur, mu_T = M_i$mu_T, loglik = logl)))
  }
  iter_i %>% cat("\nNumber of steps:", ., "\n")
  logl <- logl[1:iter_i,]
  logl2 <- logl2[1:iter_i]
  any(diff(logl2[-1])< 0) %>% cat("Any negative increment?", ., "\n")
  tail(logl2,1) %>% cat("Log likelihood:", ., "\n")
  tail(diff(logl2),1) %>% cat("Last increment:", ., "\n")

  #cbind(prm.cur$W, prm.tst$W) %>% print
  cat("\nSigmaT are ", svd(crossprod(prm.cur$W))$d, "\n")
  cat("Betas are ", prm.cur$beta, "\n")

  return(list(param = prm.cur, mu_T = M_i$mu_T, loglik = logl, loglik_sum = logl2))
}
