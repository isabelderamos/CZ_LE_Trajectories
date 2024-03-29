# LT function from Camarda 
# https://sites.google.com/site/carlogiovannicamarda/r-stuff/life-expectancy-confidence-interval

## function for constructing a classic (& rather general) lifetable
lifetable <- function(x, Nx, Dx, sex="M", ax=NULL){
  m <- length(x)
  mx  <- Dx/Nx
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="W"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      if(sex=="B"){
        if(mx[1]>=0.107){
          ax[1] <- 0.340
        }else{
          ax[1] <- 0.049 + 2.742*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, Nx, Dx, mx, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}

## function for constructing a lifetable starting from probabilities
lifetable.qx <- function(x, qx, sex="M", ax=NULL, last.ax=5.5){
  m <- length(x)
  n <- c(diff(x), NA)
  qx[is.na(qx)] <- 0
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- last.ax
    }else{    
      if(sex=="W"){
        if(qx[1]>=0.1){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.05 + 3*qx[1]
        }
      }
      if(sex=="M"){
        if(qx[1]>=0.1){
          ax[1] <- 0.33
        }else{
          ax[1] <- 0.0425 + 2.875*qx[1]
        }
      }
      if(sex=="B"){
        if(qx[1]>=0.1){
          ax[1] <- 0.340
        }else{
          ax[1] <- 0.04625 + 2.9375*qx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- last.ax
    }
  }
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]*last.ax
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return.df <- data.frame(x, n, ax, qx, px, lx, dx, Lx, Tx, ex)
  return(return.df)
}




## x=agesA
## Nx=myE
## Dx=myD
## ax=rep(0.5,mA)
## ns=100
## level=0.90
## which.x = 31
## sex = "M"
## last.ax=0.5
## general function for getting CI of life expectancy at x
CIex <- function(x, Nx, Dx, sex = "M", ax = NULL,
                 which.x=0, ns=1000, level=0.95){
  ## point-estimated lifetable
  LT <- lifetable(x, Nx, Dx, sex = sex, ax = ax)
  ## number of ages
  m <- nrow(LT)
  ## estimated probs
  qx <- LT$qx
  ## trials for binomial, rounded
  Ntil <- round(Dx/qx)
  ## ax for last age
  last.ax <- LT$ax[m]
  ## simulated death counts
  ## from Binomial distribution
  Y <- suppressWarnings(matrix(rbinom(m*ns,
                                      Ntil,
                                      qx),
                               m,ns))
  ## simulated probabilities
  QX <- Y/Ntil
  ## which age?
  wh <- which(x==which.x)
  ## for all replicates, compute life expectancy
  ## ## by a for-loop
  ## exsimA <- rep(0,ns)
  ## for(s in 1:ns){
  ##   exsimA[s] <-lifetable.qx(x, qx=QX[,s], sex,
  ##                        last.ax=last.ax)$ex[wh]
  ## }
  ## by apply command
  ## (slighly faster, less intuitive)
  fun.ex <- function(qx){
    return(lifetable.qx(x=x, qx, sex=sex,
                        last.ax=last.ax)$ex[wh])
  }
  exsim <- apply(QX, 2, fun.ex)

  ## confidence interval
  CI <- quantile(exsim,
                 probs = c((1-level)/2,
                           1 - (1-level)/2))
  ## output
  out <- list(ex=LT$ex[wh],
              meanex=mean(exsim),
              CIex=CI,
              exsim=exsim,
              which.x=which.x)
  return(out)
}
