mvregerrGS <- function(Z,
                       Z.errV,
                       X,
                       nchain=1e2,init=list(nu=5),
                       every=10,XP=NULL,record=c('Y')){
  n <- nrow(Z); k <- ncol(Z); p <- ncol(X)
  print(c(n,p,k))
  if(is.null(init$Psi)){
    Psi <- diag(k)
  } else {Psi <- init$Psi}
  if(is.null(init$B)){
    B <- matrix(5,nrow=p,ncol=k)
  } else {B <- init$B}
  nu <- init$nu
  ## containers for keeping data
  ans <- list()
  if('Y' %in% record) ans$Y <- list()
  if('beta' %in% record) ans$beta <- list()
  if('Sig' %in% record) ans$Sig <- list()
  if(!is.null(XP)) ans$YP <- list()         #for predictions
  XtX <- t(X) %*% X
  Sig <- Psi
  iB <- 1/B
  viB <- c(iB)
  Y <- Z
  bet <- solve(diag(k) %x% XtX,c(t(X) %*% Y))
  bet <- matrix(bet,nrow=p,ncol=k)      #MLE
  for(i in 1:nchain){
    if(!i%%every)cat('i=',i,'\n')
    ## y update
    ## TODO think on doing as single matrix calc; maybe faster; my N is not big
    iSig <- solve(Sig)
    for(r in 1:n){
      is2 <- diag(1/Z.errV[r,]) + iSig
      s2 <- solve(is2)                  #TODO consider efficiency: chol for lots
      v <- s2 %*% (diag(1/Z.errV[r,]) %*% (Z[r,]) + iSig %*% t(X[r,] %*% bet))
      Y[r,] <- mvrnorm(n=1,mu=v,Sigma=s2)
    }
    ## beta update
    iXi <- iSig %x% XtX   + diag(viB)
    ## ## slow but clear version
    ## Xi <- solve(iXi)
    ## b <- Xi %*% c(t(X) %*% Y %*% iSig) #1
    ## bet <- mvrnorm(n=1,mu=b,Sigma=Xi) #single draw for vector beta
    ## bet <- matrix(bet,nrow = p,ncol=k)
    ## ## more efficient version
    uiXi <- chol(iXi)
    b <- backsolve(upper.tri=TRUE,transpose=FALSE,r=uiXi,
                   x=backsolve(upper.tri=TRUE,transpose=TRUE,r=uiXi,
                               x=c(t(X) %*% Y %*% iSig))) #see 1 above
    bet <- backsolve(upper.tri=TRUE,transpose=FALSE,
                   r=uiXi,x= matrix(rnorm(length(b)),nrow=length(b),ncol=1))
    bet <- c(bet) + c(b)
    bet <- matrix(bet,nrow = p,ncol=k)
    ## update: \Sigma | Y \beta ~ IW(\Psi + (Y-X\beta)^T(Y-X\beta),\nu + n)
    V <- Y - X %*% bet
    V <- t(V) %*% V
    Sigma <- riwish(nu + n, Psi + V) #
    ## predictions
    if(!is.null(XP)){
      ans$YP[[i]] <- t(apply(XP %*% bet,1,function(x) mvrnorm(n=1,mu=x,Sigma=Sig)))
    }
    ## store
    if('Y' %in% record) ans$Y[[i]] <- Y
    if('beta' %in% record) ans$beta[[i]] <- bet
    if('Sig' %in% record) ans$Sig[[i]] <- Sig
  }
  ans
}