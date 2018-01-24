# gridding by a linear hierarchical model
# using sparse linear solver


# this uses sparse Cholesky
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
  HH <- rbind(H1,H2,H3,H4)
  yy <- rbind(y,matrix(0.0,nrow=nrow(H2)+nrow(H3)),matrix(c(priormodel),nrow=m*n))
  R<-chol(crossprod(HH))
  x<-backsolve(R,backsolve(R,crossprod(HH,yy),transpose=TRUE),transpose=FALSE)
#  out <- modelgrid
#  out@data[,variable] <- x
  dim(x)<-c(n,m)
  return(x)
}

# Using sparse solve
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
  x <- solve(crossprod(H1) + crossprod(H2) + crossprod(H3) + crossprod(H4), crossprod(H1,y) + crossprod(H4,matrix(c(priormodel),nrow=m*n)))
  dim(x)<-c(n,m)
  x<-as.matrix(x)
  #  out <- modelgrid
  #  out@data[,variable] <- x
  return(x)
}

