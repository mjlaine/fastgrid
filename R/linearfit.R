# gridding by a linear hierarchical model
# using sparse linear solver


# this uses sparse Cholesky

#' @export 
linearinterpolationtogrid <- function(obspoints, modelgrid, sigobs=1.0, sigmodel=1.0, lam=c(1,1), d=2,
                   priormodel = 0.0, method = 'bilinear', variable = 'temperature') {
  elon <- gridlon(modelgrid)
  elat <- gridlat(modelgrid)
  m <- length(elon)
  n <- length(elat)
  londx <- abs(elon[2]-elon[1])
  latdx <- abs(elat[2]-elat[1])
  if (is(modelgrid,'SpatialGridDataFrame')) {
    # priormodel <- c(modelgrid@data[,variable])
  }
  H1 <- Hmat(elon,elat,sp::coordinates(obspoints),method=method) / sigobs
  y <- as.matrix(obspoints@data[,variable])/sigobs
  # smoothing matrices for rows and columns
  H2 <- Diagonal(n)%x%diff(Diagonal(m),differences=d) * lam[1]/londx**d
  H3 <- diff(Diagonal(n),differences=d)%x%Diagonal(m) * lam[2]/latdx**d
  H4 <- Diagonal(n*m)/sigmodel
  HH <- rbind(H1,H2,H3,H4)
  yy <- rbind(y,matrix(0.0,nrow=nrow(H2)+nrow(H3)),matrix(c(priormodel),nrow=m*n))
  R<-chol(t(HH)%*%HH)
  x<-backsolve(R,backsolve(R,t(HH)%*%yy,transpose=TRUE),transpose=FALSE)
#  out <- modelgrid
#  out@data[,variable] <- x
  dim(x)<-c(n,m)
  return(x)
}

# this is as slow or slower
#' @export 
linearinterpolationtogrid2 <- function(obspoints, modelgrid, sigobs=1.0, sigmodel=1.0, lam=c(1,1), d=2,
                                      priormodel = 0.0, method = 'bilinear', variable = 'temperature') {
  elon <- gridlon(modelgrid)
  elat <- gridlat(modelgrid)
  m <- length(elon)
  n <- length(elat)
  londx <- abs(elon[2]-elon[1])
  latdx <- abs(elat[2]-elat[1])
  if (is(modelgrid,'SpatialGridDataFrame')) {
    # priormodel <- c(modelgrid@data[,variable])
  }
  H1 <- Hmat(elon,elat,sp::coordinates(obspoints),method=method) / sigobs
  y <- as.matrix(obspoints@data[,variable])/sigobs
  # smoothing matrices for rows and columns
  H2 <- Diagonal(n)%x%diff(Diagonal(m),differences=d) * lam[1]/londx**d
  H3 <- diff(Diagonal(n),differences=d)%x%Diagonal(m) * lam[2]/latdx**d
  H4 <- Diagonal(n*m)/sigmodel
# HH <- rbind(H1,H2,H3,H4)
#  yy <- rbind(y,matrix(0.0,nrow=nrow(H2)+nrow(H3)),matrix(c(priormodel),nrow=m*n))
#  R<-chol(t(HH)%*%HH)
  R <- chol(t(H1)%*%H1 + t(H2)%*%H2 + t(H3)%*%H3 + t(H4)%*%H4)
#  x<-backsolve(R,backsolve(R,t(HH)%*%yy,transpose=TRUE),transpose=FALSE)
  x<-backsolve(R,backsolve(R,t(H1)%*%y + t(H4)%*%matrix(c(priormodel),nrow=m*n),transpose=TRUE),transpose=FALSE)
  #  out <- modelgrid
  #  out@data[,variable] <- x
  dim(x)<-c(n,m)
  return(x)
}


