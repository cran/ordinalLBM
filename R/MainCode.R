library(RColorBrewer)


# ********************************************************************************************************
# The next lines include some functions for internal use only. You do never need to call them explicitly.

## -- M-EM for the block pair (q,r), not involving the LBM part 
ExpMax <- function(y, K, gammas, lgammas, R, eps, it_max, mu = 0, s2 = 1){
  N <- length(y)
  
  s <- sqrt(s2)
  old_lik = -Inf
  new_lik = -1
  
  it_count <- 0
  
  while (abs(new_lik - old_lik)>eps && it_count < it_max){
    it_count <- it_count + 1
    old_lik = new_lik;
    
    ## E-step
    tmp_1 <- pnorm((gammas - mu)/s)    # pnorm(betas)
    tmp_2 <- pnorm((lgammas - mu)/s)   # pnorm(alphas)
    
    eta <- tmp_1 - tmp_2
    eta[ eta<.Machine$double.xmin ] <- .Machine$double.xmin  # no numerical issues!
    
    # 
    tmp_1 <- dnorm((gammas - mu)/s)   # dnorm(betas)
    tmp_2 <- dnorm((lgammas - mu)/s)  # dnorm(alphas)
    
    deta <- tmp_1 - tmp_2
    
    # self made calculations 
    m <- mu - sqrt(s2)*(deta/eta)
    m2 <- vector("numeric", length = K)
    for (idx in 1:K){
      beta = min((gammas[idx] - mu)/s, .Machine$double.xmax)
      alpha = max((lgammas[idx] - mu)/s, -.Machine$double.xmax)
      m2[idx] = mu^2 - 2*s*mu*(deta[idx]/eta[idx]) - s2*(beta*dnorm(beta) - alpha*dnorm(alpha))/eta[idx] + s2
    }
    
    ## M-step
    Yf <- y[y!=0]
    Nf <- length(Yf)
    tmp_3 <- table(Yf)
    if (length(tmp_3)!=K){
      nms <- as.numeric(names(tmp_3))
      sq <- c(1:K)
      pos <- sq %in% nms
      pois <- rep(0, times = K)
      pois[pos] <- as.numeric(table(Yf))
      pois <- pois/Nf
    }
    else pois <- as.numeric(tmp_3)/Nf
    mu <- sum(m*pois)
    s2 <- sum(m2*pois) - mu^2
    s <- sqrt(s2)  
    
    new_lik <- LikG(y, eta)
  }
  return(list(mu, s2, eta))
}
## -- Greedy Classification step
Greedy <- function(Y, PI, eta_mat, R, C, rho, delta){
  nb_row_swaps <- 0
  nb_col_swaps <- 0
  N <- nrow(Y)
  P <- ncol(Y)
  Q <- dim(eta_mat)[1]
  L <- dim(eta_mat)[2]
  K <- dim(eta_mat)[3]
  ### *** row swaps ***
  for (idx in 1:N){
    #print(idx)
    # test swap
    gain <- 0
    best <- R[idx]
    if (sum(R==best)>1){
      for (q in 1:Q){
        tmp_gain <- log(rho[q]) - log(rho[R[idx]])
        for (idy in 1:P){
          if(Y[idx, idy]!=0) tmp_gain <- tmp_gain + log(eta_mat[q, C[idy], (Y[idx, idy])]) - log(eta_mat[R[idx], C[idy], (Y[idx, idy])]) + log(PI[q, C[idy]]/PI[R[idx], C[idy]])
          else tmp_gain <- tmp_gain + log((1 - PI[q,C[idy]])/(1 - PI[R[idx], C[idy]]))
        }
        if (tmp_gain>gain){
          gain <- tmp_gain
          best <- q         
        }
      }
      # do swap
      if (best != R[idx]){
        R[idx] = best
        nb_row_swaps <- nb_row_swaps + 1
      }
    }
  }
  ### *** col swaps ***
  for (idx in 1:P){
    # test swap
    gain <- 0
    best <- C[idx]
    if (sum(C==best)>1){
      for (l in 1:L){
        tmp_gain <- log(delta[l]) - log(delta[C[idx]])
        for (idy in 1:N){
          if(Y[idy, idx]!=0) tmp_gain <- tmp_gain + log(eta_mat[R[idy], l, (Y[idy,idx])]) - log(eta_mat[R[idy], C[idx], (Y[idy,idx])]) + log(PI[R[idy], l]/PI[R[idy], C[idx]])
          else tmp_gain <- tmp_gain + log((1 -  PI[R[idy], l])/(1 - PI[R[idy], C[idx]]))
        }
        if (tmp_gain>gain){
          gain <- tmp_gain
          best <- l
        }
      }
      # do swap
      if (best != C[idx]){
        C[idx] = best
        nb_col_swaps <- nb_col_swaps + 1
      }
    }
  }
  return(list(R, C, nb_row_swaps, nb_col_swaps))
}
## -- To simulate notes
GetNote <- function(v, x){
  i <- 1
  while (x > v[i]){
    i = i+1
    if (i > length(v)) break; 
  }
  return(i-1)
}
## -- computing the log-likelihood (binary LBM) for each block pair (q,r)
LikB <- function(piqr, ActualInt, PotentialInt){
  return(ActualInt*log(piqr/(1 - piqr)) + PotentialInt*log(1 - piqr))
}
## -- computing the log-likelihood (latent Gaussian random variables) for each block pair (q,r)
LikG <- function (yE, eta){
  if (sum(yE!=0) == 0){
    warning(" Warning: All the scores in the block are equal to zero! ")
  }
  y <- yE[yE!=0]       ## removing missing data (if any)
  K <- length(eta)
  y_ <- table(y)
  if (length(y_) !=0){
    nms <- as.numeric(names(y_))
    sq <- c(1:K)
    pos <- sq %in% nms
    tmp <- rep(0, times = K)
    tmp[pos] <- as.numeric(y_)
    y_ <- tmp
  }
  else y_ <- as.numeric(y_)
  return(sum(y_*log(eta)))
}
### Global log-likelihood
GLik <- function(Y, R, C, PI, eta_mat, rho, delta){
  out <- 0
  for (idx in 1:nrow(Y)){
    for (idy in 1:ncol(Y)){
      if(Y[idx,idy]!=0) out = out + log(eta_mat[R[idx], C[idy], (Y[idx, idy])]) + log(PI[R[idx], C[idy]]) 
      else out = out + log(1-PI[R[idx], C[idy]]) 
    }
  }
  for (idx in 1:nrow(Y)) out = out + log(rho[R[idx]]) 
  for (idy in 1:ncol(Y)) out = out + log(delta[C[idy]])
  return(out)
}
### To avoid degenerated cases with empty classes
PiCheck <- function(piqr){
  if (piqr < .Machine$double.eps) return(.Machine$double.eps)
  else if (piqr == 1) return(1-.Machine$double.eps)
  else return(piqr)
}

# *********************************************************************************************************
# *********************************************************************************************************
# ---  Main functions  to simulate, fit to data and plot OLBM  

# 2.

#' Fitting OLBM to the data
#' 
#' It estimates the OLBM model parameters as well as the most likely posterior cluster assignments by maximum likelihood.
#' @param Y An M x P ordinal matrix, containing ordinal entries from 1 to K. Missing data are coded as zeros.
#' @param Q The number of row clusters.
#' @param L The number of column clusters. 
#' @param init A string specifying the initialisation type. It can be "kmeans" (the default) or "random" for a single random initialisation.
#' @param eps When the difference between two consecutive vaules of the log-likelihood is smaller than eps, the M-EM algorithms will stop.
#' @param it_max The maximum number of iterations that the M-EM algorithms will perform (although the minimum tolerance eps is not reached). 
#' @param verbose A boolean specifying whether extended information should be displayed or not (TRUE by default).
#' 
#' @return It returns an S3 object of class "olbm" containing 
#'           \item{estR}{the estimated row cluster memberships.}
#'           \item{estC}{the estimated column cluster memberships.}
#'           \item{likeli}{the final value of the log-likelihood.}
#'           \item{icl}{the value of the ICL criterion.}
#'           \item{Pi}{the Q x L estimated connectivity matrix.}
#'           \item{mu}{a Q x L matrix containing the estimated means of the latent Gaussian distributions.}
#'           \item{sd}{a Q x L matrix containing the estimated standard deviations of the latent Gaussian distributions.}
#'           \item{eta}{a Q x L x K array whose entry (q,l,k) is the estimated probability that one user in the q-th row cluster assign the score k to one product in the l-th column cluster.}
#'           \item{rho}{the estimated row cluster proportions.}
#'           \item{delta}{the estimated column cluster proportions.}
#'           \item{initR}{the initial row cluster assignments provided to the C-EM algorithm.}
#'           \item{initC}{the initial column cluter assignments provided to the C-EM algorigthm.}
#'           \item{Y}{the input ordinal matrix Y.}
#'           \item{thresholds}{the values (1.5, 2.5, ... , K-0.5) of the thresholds, defined inside the function olbm.}
#' @export 
#' @references
#' Corneli M.,Bouveyron C. and Latouche P. (2019) \emph{Co-Clustering of ordinal data via latent continuous random variables and a classification EM algorithm. }(\url{https://hal.archives-ouvertes.fr/hal-01978174})
#' @examples
#' 
#' data(olbm_dat)
#' res <- olbm(olbm_dat$Y, Q=3, L=2)                       
olbm <- function(Y, Q, L, 
                 init = "kmeans",  
                 eps = 10e-5, 
                 it_max =  500,
                 verbose = TRUE){
  
  K <- max(Y)
  N <- nrow(Y)
  P <- ncol(Y)
  
  if (init == "random"){
    # -- Random init
    R <- sample(seq(1:Q), N, replace = TRUE, prob = rep(1, times = Q)) 
    C <- sample(seq(1:L), P, replace = TRUE, prob = rep(1, times = L))
    initR <- R
    initC <- C
  }
  
  else if (init == "kmeans"){
    # -- Kmeans init
    R <- kmeans(Y, Q)
    R <- R$cluster
    C <- kmeans(t(Y), L)
    C <- C$cluster
    initR <- R
    initC <- C
  }
  
  else{
    warning(" Wrong assignment to the variable init! ")
  }
  
  Y_ <- Y[Y!=0]
  
  # thresholds
  gammas <- 0.5 + c(1:(K-1))
  gammas <- c(gammas,Inf)
  lgammas <- c(-Inf, gammas)
  lgammas <- lgammas[-length(lgammas)]
  
  nb_row_swaps <- nb_col_swaps <- -1       # at random, just to be different from zero!
  # mu <- matrix(0, nrow = Q, ncol = L)
  # s2 <- matrix(1, nrow = Q, ncol = L)
  mu <- matrix(mean(lgammas[-1]), nrow = Q, ncol = L)
  if ( length(lgammas) < 3 ) s2 <- matrix(1, nrow = Q, ncol = L)
  else s2 <- matrix(sd(lgammas[-1]), nrow = Q, ncol = L)
  eta_mat <- array(dim = c(Q,L,K))
  PI <- matrix(0, nrow = Q, ncol = L) 
  while(nb_row_swaps!=0 || nb_col_swaps!=0){
    if (verbose) print("C-EM: ")
    ## *** best rho
    tmp <- table(R)
    if ( length(tmp)!=Q){
      nms <- as.numeric(names(tmp))
      sq <- c(1:Q)
      pos <- sq %in% nms
      tmp_ <- rep(0, times = Q)
      tmp_[pos] <- as.numeric(tmp)
      tmp <- tmp_
    }
    rho <- as.numeric(tmp)/length(R)
    rho[rho == 0] <- .Machine$double.xmin
    ## *** best delta
    tmp <- table(C)
    if ( length(tmp)!=L){
      nms <- as.numeric(names(tmp))
      sq <- c(1:L)
      pos <- sq %in% nms
      tmp_ <- rep(0, times = L)
      tmp_[pos] <- as.numeric(tmp)
      tmp <- tmp_
    }
    delta <- as.numeric(tmp)/length(C)
    delta[delta == 0] <- .Machine$double.xmin
    
    global_lik <- sum(log(rho[R])) + sum(log(delta[C]))
    for (q in 1:Q){
      for (l in 1:L){
        yql_ <- as.numeric(Y[R==q, C==l])
        yql <- yql_[yql_ != 0]
        if (length(yql)!=0){ # In the block pair (q,l) there are some interactions 
          ActualInt <- length(yql)
          PotentialInt <- length(yql_)
          PI[q,l] <- ActualInt/PotentialInt
          PI[q,l] <- PiCheck(PI[q,l])
          out <- ExpMax(yql, K, gammas, lgammas, R, mu = mu[q,l], s2 = s2[q,l], eps = eps, it_max = it_max)
          mu[q,l] <- out[[1]]
          s2[q,l] <- out[[2]]
          eta_mat[q,l,] <- out[[3]]
          #if (sum(is.nan(eta_mat[q,l,])) != 0) eta_mat[q,l,] <- rep(1/K, times = K)
          global_lik <- global_lik + LikB(PI[q,l], ActualInt, PotentialInt) + LikG(yql,eta_mat[q,l,])
        }
        else{ # there are no interactions in the block pair (q,l)
          ActualInt <- 0
          PotentialInt <- length(yql_)
          PI[q,l] <- ActualInt/PotentialInt
          PI[q,l] <- PiCheck(PI[q,l])
          global_lik <- global_lik + LikB(PI[q,l], ActualInt, PotentialInt) #+ LikG(yql,eta_mat[q,l,])
        }
      }
    }
    if (verbose){
      cat("Global likelihood: ", global_lik, "\n")
      print("Greedy: ")
    }
    out <- Greedy(Y, PI, eta_mat, R, C, rho, delta)
    R <- out[[1]]
    C <- out[[2]]
    nb_row_swaps <- out[[3]]
    nb_col_swaps <- out[[4]]
  }
  ## Last M-EM
  global_lik <- sum(log(rho[R])) + sum(log(delta[C]))
  for (q in 1:Q){
    for (l in 1:L){
      yql_ <- as.numeric(Y[R==q, C==l])
      yql <- yql_[yql_ != 0]
      if (length(yql)!=0){ # In the block pair (q,l) there are some interactions 
        ActualInt <- length(yql)
        PotentialInt <- length(yql_)
        PI[q,l] <- ActualInt/PotentialInt
        PI[q,l] <- PiCheck(PI[q,l])
        out <- ExpMax(yql, K, gammas, lgammas, R, mu = mu[q,l], s2 = s2[q,l], eps = eps, it_max = it_max)
        mu[q,l] <- out[[1]]
        s2[q,l] <- out[[2]]
        eta_mat[q,l,] <- out[[3]]
        #if (sum(is.nan(eta_mat[q,l,])) != 0) eta_mat[q,l,] <- rep(1/K, times = K)
        global_lik <- global_lik + LikB(PI[q,l], ActualInt, PotentialInt) + LikG(yql,eta_mat[q,l,])
      }
      else{ # there are no interactions in the block pair (q,l)
        ActualInt <- 0
        PotentialInt <- length(yql_)
        PI[q,l] <- ActualInt/PotentialInt
        PI[q,l] <- PiCheck(PI[q,l])
        global_lik <- global_lik + LikB(PI[q,l], ActualInt, PotentialInt) #+ LikG(yql,eta_mat[q,l,])
      }
    }
  }
  
  ## ICL criterion
  icl <- global_lik - Q*L*log(length(Y_)) - Q*L*log(N*P)/2 - Q*log(N)/2 - L*log(P)/2
  
  out_list <- list(estR=R, estC=C, likeli=global_lik, icl=icl, Pi=PI, mu=mu,
                   Sigma = s2, eta = eta_mat, rho = rho, delta = delta, initR = initR, 
                   initC=initC, Y = Y, thresholds = gammas[-length(gammas)])
  class(out_list) <- "olbm"
  return(out_list)
}

# 1.

#' Simulate OLBM data
#' 
#' It simulates an ordinal data matrix according to OLBM.
#' @param M The number of rows of the ordinal matrix Y.
#' @param P The number of columns of the ordinal matrix Y.
#' @param Pi A Q x L connectivity matrix to manage missing data (coded az zeros in Y).
#' @param rho  A vector of length Q, containing multinomial probabilities for row cluster assignments.
#' @param delta A vector of length L, containing multinomial probabilities for column cluster assignments.
#' @param mu A Q x L matrix containing the means of the latent Gaussian distributions.
#' @param sd A Q x L matrix containing the standard deviations of the latent Gaussian distributions.
#' @param thresh A K+1 vector containing the sorted tresholds used to simulate the ordinal entries in Y, where K is the number of ordinal modalities. The first entry in tresh must be -Inf, the last entry +Inf.
#'
#' @return   It returns a list containing:  
#'               \item{Y}{An M x P matrix. The observed ordinal entries are integers between 1 and K. Missing data are coded as zeros.}
#'               \item{Rclus}{A vector of length M containing the row cluster memberships.}
#'               \item{Cclus}{A vector of length P containing the column cluster memberships.}
#' @export     
#' @references 
#' Corneli M.,Bouveyron C. and Latouche P. (2019) \emph{Co-Clustering of ordinal data via latent continuous random variables and a classification EM algorithm. }(\url{https://hal.archives-ouvertes.fr/hal-01978174})
#' @examples
#' 
#' M <- 150                                    
#' P <- 100 
#' Q <- 3
#' L <- 2
#' 
#' ## connectivity matrix
#' Pi <- matrix(.7, nrow = Q, ncol = L)
#' Pi[1,1] <- Pi[2,2] <- Pi[3,2] <- .5
#' 
#' ## cluster memberships proportions
#' rho <- c(1/3, 1/3 ,1/3)
#' delta <- c(1/2, 1/2)
#' 
#' # Thresholds
#' thresh <- c(-Inf, 2.37, 2.67, 3.18, 4.33, Inf)     # K = 5
#'
#' ## Gaussian parameters
#' mu <- matrix(c(0, 3.4, 2.6, 0, 2.6, 3.4), nrow = Q, ncol = L)   
#' sd <- matrix(c(1.2,1.4,1.0,1.2,1.4,1.0), nrow = Q, ncol = L)
#' 
#' ## Data simulation
#' dat <- simu.olbm(M, P, Pi, rho, delta, mu, sd, thresh)
#' 
simu.olbm <- function(M,P,Pi,rho,delta,mu,sd,thresh){
  Q <- nrow(Pi)
  L <- ncol(Pi)
  Rclus <- sample(c(1:Q), size = M, replace = TRUE, prob = rho)
  Cclus <- sample(c(1:L), size = P, replace = TRUE, prob = delta)
  ## simulating the latent vbs and the notes
  Y <- Z <- matrix(0,nrow = M, ncol = P)
  for (i in 1:M){
    for (j in 1:P){
      u <- runif(1)
      if (u < Pi[Rclus[i], Cclus[j]]){
        q = Rclus[i]
        r = Cclus[j]
        Z[i,j] <- rnorm(1, mean = mu[q,r], sd = sd[q,r])
        Y[i,j] <- GetNote(thresh, Z[i,j])
      }
    }
  }
  return(list(Y=Y, Rclus = Rclus, Cclus = Cclus))
}
 
#3.

#' Plot OLBM
#' 
#' It plots the re-organized incidence matrix and/or the estimated Gussian densities.
#' @param x    The "olbm" object output of the function olbm.
#' @param type A string specifying the type of plot to be produced. The currently supported values are "hist" and "incidence".
#' @param ...  Additional parameters to pass to sub-functions.
#' @export
#' @import
#' grDevices
#' graphics
#' stats
#'  
#' @examples
#' data(olbm_dat)
#' res <- olbm(olbm_dat$Y, Q=3, L=2)   
#' plot(res, "hist")
#' plot(res, "incidence")
plot.olbm <- function(x, type = "hist", ...){
  res <- x
  paletta <- RColorBrewer::brewer.pal(12, name = "Paired")
  palette(paletta)
  Rclus <- res$estR
  Cclus <- res$estC
  Y <- res$Y
  M <- nrow(Y)
  P <- ncol(Y)
  K <- max(Y)
  Q <- max(Rclus)
  L <- max(Cclus)
  mu <- res$mu
  s2 <- res$Sigma
  delta <- res$delta
  rho <- res$rho
  Pi <- res$Pi
  if (type == "hist"){
    #op <- par(mfrow = c(Q,L),mar=c(2,2,2,2),oma=c(0,1,1,1), bty = "n")
    griglia <- seq(from = 0, to = (K+1), by = 0.01)
    idx <- 0
    op <- par(mfrow = c(Q,L),mar=c(2,2,2,2),oma=c(0,1,1,1), bty = "n")
    for (q in 1:Q){
      for (l in 1:L){
        t1 <- sum(Rclus == q)
        t2 <- sum(Cclus == l)
        selection <- Y[Rclus==q, Cclus==l]
        selection <- selection[selection!=0]
        idx <- idx + 1
        hist(selection, probability = TRUE, main = "", col = idx, 
             ylim = c(0,1.2), breaks = seq(0.5, (K+0.5), by =1), ...)
        mn <- 100 * t1*t2/(M*P)
        if (sum(selection) != 0){
          valori <- dnorm(griglia, mean = mu[q,l], sd = sqrt(s2[q,l]))
          points(griglia, valori, type = "l", col = "darkgrey", lwd = 2, ...) #, ylim = c(0,0.9))
        }
        qc <- as.character(q)
        lc <- as.character(l)
        qclc <- paste(q,l,sep ="")
        mtext(bquote(rho[.(qc)]*delta[.(lc)] == .(round(100*rho[q]*delta[l], digits = 2)) ~'%'*"    "*pi[.(qclc)] == .(round(100*Pi[q,l], digits = 2)) ~'%'), cex = 0.65)
      }
    }
  }
  if (type == "incidence"){
    op <- par(mfrow = c(1,1), bty = "n")
    I <- Y
    I[I != 0] <- 1
    sR <- seq(1:M)
    Rclus.copy <- Rclus
    names(Rclus.copy) <- sR
    Rclus.copy <- sort(Rclus.copy)
    posR <- as.numeric(names(Rclus.copy))
    sC <- seq(1:P)
    Cclus.copy <- Cclus
    names(Cclus.copy) <- sC
    Cclus.copy <- sort(Cclus.copy)
    posC <- as.numeric(names(Cclus.copy))
    J <- I[posR, posC]
    image(t(J)[,M:1], col=gray(255:0/255), useRaster = TRUE, xaxt = 'n', yaxt = 'n', ...)  
    # hlines
    lev_hl <- cumsum(rev(table(Rclus.copy)))/M
    lev_hl <- lev_hl[-Q]
    abline(h = lev_hl, col = "red", lwd = 1.5)
    # print the number of  row cluster
    lag_lev_hl <- c(0, lev_hl)
    lev_hl <- c(lev_hl,1)
    pos <- 0.5*(lev_hl - lag_lev_hl)
    pos <- lag_lev_hl + pos
    axis(side = 2, at = pos, labels = c(Q:1), tick = FALSE)
    # vlines
    lev_vl <- cumsum(table(Cclus.copy))/P
    lev_vl <- lev_vl[-L]
    abline(v = lev_vl, col = "red", lwd = 1.5)
    # print the number of  col cluster
    lag_lev_vl <- c(0, lev_vl)
    lev_vl <- c(lev_vl,1)
    pos <- 0.5*(lev_vl - lag_lev_vl)
    pos <- lag_lev_vl + pos
    axis(side = 3, at = pos, labels = c(1:L), tick = FALSE)
  }
}


