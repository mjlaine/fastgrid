### Fast Krieging by calling fortran Krieging code from R
### This only works with MOS temperature data in predefined format
###
### compilation of the shared code by a command like:
### R CMD SHLIB kriegecode.f90 -lRblas -lRlapack
###
### or install as package:
### R CMD INSTALL fastgrid_n.n.tar.gz

### marko.laine@fmi.fi

## TODO: error processing!!!

#' fastgrid: faster Kriging by some f90 code
#'
#' Thre are some functions here to help with gridding point data to a regular grid.
#' 
#' @section functions:
#' See \code{\link{fastkriege}}.
#'
#' @docType package
#' @name fastgrid
NULL

# methods package needed if run by Rscript
# but seems to need library(methods) in the file.R too (for older version of Matrix)
#' @import sp Matrix
NULL

## This loads the code to R memory
# if (!is.loaded("kriegecode")) {
#    dyn.load("kriegecode.so")
#}

#require(Matrix) # for sparseMatrix
#require(sp) # for spatial objects

## to replace gstat kriege command for prediction over a grid
#' Calculate Kriging predictions
#' 
#' @param trend_model formula to build the regression matrices
#' @param data data containing station temperatures as SpatialPoints
#' @param grid definition of model grid, must contain longitude, latitude coordinates and variables used in trend_model
#' @param cov.pars c(sigmasq,phi,tausq)
#' @param bg optional background field of the same dimensions as the grid
#' @param lsm,lsmy land sea masks for grid and data
#' @param alt,alty altitude information for the grid and data
#' @param altlen range parameter for altitude (meters)
#' @param variable the name of the variable to be gridded, default "temperature"
#' 
#' @description Now assumes that the coordinate names are \code{longitude} and \code{latitude}.
#'
#' @export
fastkriege <- function(trend_model = temperature ~ -1, data, grid, cov.pars, 
                       lsm=NULL,lsmy=NULL, alt=NULL, alty=NULL, altlen=200.0,
                       bg=NULL, variable="temperature" ) {
  ## build input matrices
  ## assumes data and bg have "longitude", "latitude", "temperature", and data and grid also trend model variables
  
  trend_model_noy <- trend_model[-2]
  grid.variables <- all.vars(trend_model_noy)
  
  if (!is.null(bg)) {
    ## map from from bg grid to station locations
    ## need the grid for obs operator, FIX me
    s<-sp::summary(grid)$grid
    elon<-seq.int(from=s[1,1],by=s[1,2],len=s[1,3])
    elat<-seq.int(from=s[2,1]+s[2,2]*(s[2,3]-1) ,by=-s[2,2],len=s[2,3])
    H<-f90Hmat(elon,elat,cbind(data$longitude,data$latitude))
    mu <- as.matrix(H%*%as.matrix(bg@data[,variable]))
  }
  else {
    mu <- 0.0
  }
  
  ## non missing of all trend model variables
  igrid <- complete.cases(as.data.frame(grid)[,grid.variables])
  
  B <- buildcovmat(coords=coordinates(data), cov.model = "exp", cov.pars=cov.pars)
  X <- model.matrix(trend_model_noy, data=data)
  y <- as.matrix(data@data[,variable]-mu)
  predgrid<-model.matrix(trend_model_noy,data=grid)
  cy <- as.matrix(coordinates(data))
  cgrid <- as.matrix(coordinates(grid))[igrid,]
  
  ## Kriging by fortran code
  t1<-proc.time()
  if (is.null(lsm)) {
    ypred<-f90kriege(X,y,B,predgrid,cy,cgrid,cov.pars)
  }
  else if (is.null(alt) | is.null(alty)) {
    if (is.null(lsmy)) {
      lsmy <- rep(1,length(y))
    }
    B <- fixseapointsincov(B,lsmy)
    ypred<-f90kriege2(X,y,B,predgrid,lsm,lsmy,cy,cgrid,cov.pars)
  }
  else {
    if (is.null(lsmy)) {
      lsmy <- rep(1,length(y))
    }
    B <- fixseapointsincov(B,lsmy)
    # altitude based covariance matrix as in f90 code !!
#    B <- B*exp(-as.matrix((dist(alty)/altlen)^2))
    B <- B*exp(-as.matrix((dist(alty)/altlen)))
    ypred<-f90kriege3(X,y,B,predgrid,lsm,lsmy,alt,alty,cy,cgrid,c(cov.pars[1],cov.pars[2],altlen))
  }
  t2<-proc.time()-t1
  
  ypredgrid<-double(length=nrow(grid))*NA
  ypredgrid[igrid]<-ypred
  
  ## add bg field to the results
  if (is.null(bg))
    ypred2<-data.frame(temperature=ypredgrid)
  else
    ypred2<-data.frame(temperature=ypredgrid+bg@data[,variable])
  ## change the name according to the ´variable´
  names(ypred2) <- variable
  
  names(ypred2) <- c(variable)
  coordinates(ypred2)<-coordinates(grid)
  proj4string(ypred2)<-proj4string(grid)
  gridded(ypred2)<-TRUE
  
  if (!is.null(bg))
    ypred2$diff <- as.vector(ypred2@data[,variable] - bg@data[,variable])
  
  attr(ypred2,'failed') <- attr(ypred,'failed')
  
  return(ypred2)
}

#' @export
fastkriege_dev <- function(trend_model = temperature ~ -1, data, grid, cov.pars, 
                       lsm=NULL,lsmy=NULL, alt=NULL, alty=NULL, altlen=200.0,
                       bg=NULL, variable="temperature" , LapseRate = 0.0) {
  ## build input matrices
  ## assumes data and bg have "longitude", "latitude", "temperature", and data and grid also trend model variables
  
  trend_model_noy <- trend_model[-2]
  grid.variables <- all.vars(trend_model_noy)
  
  if (!is.null(bg)) {
    ## map from from bg grid to station locations
    ## need the grid for obs operator, FIX me
#    s<-sp::summary(grid)$grid
#    elon<-seq.int(from=s[1,1],by=s[1,2],len=s[1,3])
#    elat<-seq.int(from=s[2,1]+s[2,2]*(s[2,3]-1) ,by=-s[2,2],len=s[2,3])
#    H<-f90Hmat(elon,elat,cbind(data$longitude,data$latitude))
#    mu <- as.matrix(H%*%as.matrix(bg@data[,variable]))
    
    mu <- grid2points_test(bg,data,modelgrid = grid, LapseRate = LapseRate)
#    mu <- grid2points_test2(bg,data)
    
  }
  else {
    mu <- 0.0
  }
  
  ## non missing of all trend model variables
  igrid <- complete.cases(as.data.frame(grid)[,grid.variables])
  
  B <- buildcovmat(coords=coordinates(data), cov.model = "exp", cov.pars=cov.pars)
  X <- model.matrix(trend_model_noy, data=data)
  y <- as.matrix(data@data[,variable]-mu)
  predgrid<-model.matrix(trend_model_noy,data=grid)
  cy <- as.matrix(coordinates(data))
  cgrid <- as.matrix(coordinates(grid))[igrid,]
  
  ## Kriging by fortran code
  t1<-proc.time()
  if (is.null(lsm)) {
    ypred<-f90kriege(X,y,B,predgrid,cy,cgrid,cov.pars)
  }
  else if (is.null(alt) | is.null(alty)) {
    if (is.null(lsmy)) {
      lsmy <- rep(1,length(y))
    }
    B <- fixseapointsincov(B,lsmy)
    ypred<-f90kriege2(X,y,B,predgrid,lsm,lsmy,cy,cgrid,cov.pars)
  }
  else {
    if (is.null(lsmy)) {
      lsmy <- rep(1,length(y))
    }
    B <- fixseapointsincov(B,lsmy)
    # altitude based covariance matrix as in f90 code !!
    #    B <- B*exp(-as.matrix((dist(alty)/altlen)^2))
    B <- B*exp(-as.matrix((dist(alty)/altlen)))
    ypred<-f90kriege3(X,y,B,predgrid,lsm,lsmy,alt,alty,cy,cgrid,c(cov.pars[1],cov.pars[2],altlen))
  }
  t2<-proc.time()-t1
  
  ypredgrid<-double(length=nrow(grid))*NA
  ypredgrid[igrid]<-ypred
  
  ## add bg field to the results
  if (is.null(bg))
    ypred2<-data.frame(temperature=ypredgrid)
  else
    ypred2<-data.frame(temperature=ypredgrid+bg@data[,variable])
  ## change the name according to the ´variable´
  names(ypred2) <- variable
  
  names(ypred2) <- c(variable)
  coordinates(ypred2)<-coordinates(grid)
  proj4string(ypred2)<-proj4string(grid)
  gridded(ypred2)<-TRUE
  
  if (!is.null(bg))
    ypred2$diff <- as.vector(ypred2@data[,variable] - bg@data[,variable])
  
  attr(ypred2,'failed') <- attr(ypred,'failed')

  attr(ypred2,'MOSvalues') <-  as.matrix(data@data[,variable])
  attr(ypred2,'mappedvalues') <-  mu
  attr(ypred2,'stationlocations') <- cy
    
  return(ypred2)
}




## cover function for fortran kriging code in kriegecode.f90
#' @useDynLib fastgrid kriegepred 
f90kriege <- function(x,y,b,grid,cy,cgrid,covpars) {
    nobs <- nrow(x)
    npar <- ncol(x)
    ngrid <- nrow(grid)

    ## There is very little checks in the Fortran code, so check the sanity of the inputs here 
    if (nrow(y) != nobs | ncol(y) != 1 | nrow(cy) != nobs
        | ncol(grid) != npar | nrow(cgrid) != ngrid
        | ncol(cy) != 2 | ncol(cgrid) != 2
        | !is.numeric(covpars) | length(covpars) < 2) {
        stop("input dimensions do not match")
    }
    ## Call the external f90 code for Kriging predictions    
    ypred <- .Fortran("kriegepred",
                      x=as.double(x),
                      y=as.double(y),
                      b=as.double(b),
                      grid=as.double(grid),
                      ypred=double(length=ngrid),
                      cy=as.double(cy),
                      cgrid=as.double(cgrid),
                      nobs=as.integer(nobs),
                      npar=as.integer(npar),
                      ngrid=as.integer(ngrid),
                      covpars=as.double(covpars)
                      )$ypred

    if (b[1,1] == -1.0) {
        stop("Problems with the covariance information, no prediction done")
    }
    
    return(ypred)
}


## cover function for fortran kriging code in kriegecode.f90
#' @useDynLib fastgrid kriegepred2
f90kriege2 <- function(x,y,b,grid,lsm,lsmy,cy,cgrid,covpars) {
  nobs <- nrow(x)
  npar <- ncol(x)
  ngrid <- nrow(grid)
  
  ## There is very little checks in the Fortran code, so check the sanity of the inputs here 
  if (nrow(y) != nobs | ncol(y) != 1 | nrow(cy) != nobs
      | ncol(grid) != npar | nrow(cgrid) != ngrid | length(lsm) != ngrid | length(lsmy) != nobs
      | ncol(cy) != 2 | ncol(cgrid) != 2
      | !is.numeric(covpars) | length(covpars) < 2) {
    stop("input dimensions do not match")
  }
  ## Call the external f90 code for Kriging predictions    
  ypred <- .Fortran("kriegepred2",
                    x=as.double(x),
                    y=as.double(y),
                    b=as.double(b),
                    grid=as.double(grid),
                    ypred=double(length=ngrid),
                    lsm=as.double(lsm),
                    lsmy=as.double(lsmy),
                    cy=as.double(cy),
                    cgrid=as.double(cgrid),
                    nobs=as.integer(nobs),
                    npar=as.integer(npar),
                    ngrid=as.integer(ngrid),
                    covpars=as.double(covpars)
  )$ypred
  
  if (b[1,1] == -1.0) {
    stop("Problems with the covariance information, no prediction done")
  }
  
  return(ypred)
}

## cover function for fortran kriging code in kriegecode.f90
#' @useDynLib fastgrid kriegepred3
f90kriege3 <- function(x,y,b,grid,lsm,lsmy,alt,alty,cy,cgrid,covpars) {
  nobs <- nrow(x)
  npar <- ncol(x)
  ngrid <- nrow(grid)
  
  ## There is very little checks in the Fortran code, so check the sanity of the inputs here 
  if (nrow(y) != nobs | ncol(y) != 1 | nrow(cy) != nobs
      | ncol(grid) != npar | nrow(cgrid) != ngrid | length(lsm) != ngrid | length(lsmy) != nobs
      | length(alt) != ngrid | length(alty) != nobs
      | ncol(cy) != 2 | ncol(cgrid) != 2
      | !is.numeric(covpars) | length(covpars) < 3) {
    stop("input dimensions do not match")
  }
  ## Call the external f90 code for Kriging predictions    
  ypred <- .Fortran("kriegepred3",
                    x=as.double(x),
                    y=as.double(y),
                    b=as.double(b),
                    grid=as.double(grid),
                    ypred=double(length=ngrid),
                    lsm=as.double(lsm),
                    lsmy=as.double(lsmy),
                    alt=as.double(alt),
                    alty=as.double(alty),
                    cy=as.double(cy),
                    cgrid=as.double(cgrid),
                    nobs=as.integer(nobs),
                    npar=as.integer(npar),
                    ngrid=as.integer(ngrid),
                    covpars=as.double(covpars)
  )$ypred
  
  # this does not work, please fix me!!!
  if (b[1,1] == -1.0) {
    attr(ypred,'failed') <- 1
    warning("Problems with the covariance information, no prediction done")
  } else {
    attr(ypred,'failed') <- 0
  }
  
  return(ypred)
}

# Hmat is used like Hmat(lon,lat) in the code! CHANGE names!!!

### cover function for fortran code for observation operator
#' @useDynLib fastgrid obsoper
#' @export
f90Hmat <- function(mlat,mlon,obs) {
  nobs <- nrow(obs)
  nlat <- length(mlat)
  nlon <- length(mlon)
  nmod <- nlat*nlon

  if (nlat < 1 | nlon < 1 | ncol(obs) != 2 | !is.numeric(obs) )
      stop("input dimensions do not match")
  
  f90h <- .Fortran("obsoper",                   
                   mlat=as.double(mlat),
                   mlon=as.double(mlon),
                   obs=as.double(obs),
                   nlat=as.integer(nlat),
                   nlon=as.integer(nlon),
                   nobs=as.integer(nobs),
                   iind=integer(length=4*nobs),
                   jind=integer(length=4*nobs),
                   h=double(length=4*nobs)
                   )

  return(Matrix::sparseMatrix(i=f90h$iind,j=f90h$jind,x=f90h$h,dims=c(nobs,nmod)))
}



### Observation operator using bi-linear interpolation and bisection
### Returns sparse matrix
### See also: f90Hmat, which uses fortran code
#' @export
Hmat <- function(mlat,mlon,obs) {
    nobs <- nrow(obs)
    nlat <- length(mlat)
    nlon <- length(mlon)
    nmod <- nlat*nlon
    ## sparse matrix
    H<-Matrix(0,nrow=nobs,ncol=nmod,sparse=TRUE)
  
    for (i in 1:nobs) {
        olat <- obs[i,1]
        olon <- obs[i,2]
        i1i2 <- lookup(mlat,olat)
        j1j2 <- lookup(mlon,olon)
        i1<-i1i2[1]
        i2<-i1i2[2]
        j1<-j1j2[1]
        j2<-j1j2[2]
        I1 <- sub2ind(c(nlat,nlon),i1,j1)
        I2 <- sub2ind(c(nlat,nlon),i2,j1)
        I3 <- sub2ind(c(nlat,nlon),i1,j2)
        I4 <- sub2ind(c(nlat,nlon),i2,j2)
        H[i,c(I1,I2,I3,I4)] <- intcoef(olat,olon,mlat[i1],mlat[i2],mlon[j1],mlon[j2])
    }
    H
}

### utility for Hmat
sub2ind<-function(nrc,irow,icol) {
    (icol-1)*nrc[1] + irow
}

### calculate interpolation coefficients for bi-linear 2D interpolation
intcoef <- function(x1,x2,x1l,x1u,x2l,x2u) {
    X <- matrix(data=c(1,0,0,0,-1,1,0,0,-1,0,1,0,1,-1,-1,1),nrow=4,ncol=4,byrow=TRUE)
    x11 <- (x1-x1l)/(x1u-x1l)
    x22 <- (x2-x2l)/(x2u-x2l)
    h <- matrix(data=c(1,x11,x22,x11*x22),nrow=1,ncol=4)%*%X
    h
}

### calculate interpolation coefficients for nearest neighbour
intcoef.nearest <- function(x1,x2,x1l,x1u,x2l,x2u) {
  # ll,ul,lu,uu
  i <- which.min(c((x1-x1l)^2+(x2-x2l)^2, (x1-x1u)^2+(x2-x2l)^2, (x1-x1l)^2+(x2-x2u)^2, (x1-x1u)^2+(x2-x2u)^2))
  h <- matrix(data=c(0.0,0.0,0.0,0.0),nrow=1,ncol=4)
  h[1,i] = 1.0
  h
}

### table lookup by bisection 
lookup<-function(x,xi) {
    n <- length(x)
    if (x[2]-x[1] > 0) {
        l <- 1
        u <- n
        while (u-l > 1) {
            m <- floor((u+l)/2)
            if (xi>x[m]) l<-m
            else  u<-m
        }
        return(c(l,u))
    } else {
        l <- n
        u <- 1
        while (l-u > 1) {
            m <- floor((u+l)/2)
            if (xi>x[m]) l<-m
            else  u<-m
        }
        return(c(u,l))
    }
}

### calculate covariance matrix given data coordinates and spatial parameters
### only exponential covariance function implemented so far
### TODO: use spDists function from sp
buildcovmat<-function (coords = NULL, cov.model = "exp",  
                    cov.pars = stop("no cov.pars argument"), dist.method = "euclidean", ...) {
    if (is.null(coords))
        stop("coords missing")
    n <- nrow(coords)
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    if (length(cov.pars)>2) 
        nugget <- cov.pars[3]
    else
        nugget <- 0.0
    tausq <- nugget
    dists.lowertri <- as.vector(dist(coords,method=dist.method))
    if (round(1e+12 * min(dists.lowertri)) == 0) 
        warning("Two or more pairs of data at coincident (or very close) locations. \nThis may cause crashes in some matrices operations.\n")
    varcov <- matrix(0, n, n)
    if (all(sigmasq < 1e-10) | all(phi < 1e-10)) {
        varcov <- diag(x = (tausq + sum(sigmasq)), n)
    }
    else {
        covvec <- exp(-(dists.lowertri/phi))*sigmasq
        varcov[lower.tri(varcov)] <- covvec
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- covvec
        diag(varcov) <- tausq + sigmasq
    }
    return(varcov)
}

# set seapoint to landpoint covariances to zero
fixseapointsincov <- function(B,s){
  if ((dim(B)[1] != dim(B)[2]) | (dim(B)[1] != length(s)))
    stop('dimension problem in the inputs')
  sm <-  !as.logical(s%*%t(s)+(!s)%*%t(!s))
  B[sm] <- 0.0
  return(B)
}


## map grid values to points (there is one in MOSfieldutils, too)
#' @export
grid2points_test <- function(grid,data,variable="temperature", LapseRate=MOSget('LapseRate'),  modelgrid=NULL, method='bilinear') {
  grid.grid <- sp::getGridTopology(grid)
  mlon<-sp::coordinatevalues(grid.grid)$longitude
  mlat<-sp::coordinatevalues(grid.grid)$latitude
  
  nlat <- length(mlat)
  nlon <- length(mlon)
  nmod <- nlat*nlon

  # needs MOSfieldutils here  
  if (is.null(modelgrid)) {
    data("KriegeData", package = MOS.options$pkg, envir = parent.frame())
    modelgrid <- KriegeData
  }

  #  grid.alt <- grid$altitude
  grid.alt <- modelgrid$elevation
  grid.y <- as.matrix(grid@data[,variable])
  data.coord <- sp::coordinates(data)
  data.alt <- data$elevation
  
  nobs <- nrow(data)
  y <- rep(0.0,nobs)
  
  for (i in 1:nobs) {
    olat <- data.coord[i,"latitude"]
    olon <- data.coord[i,"longitude"]
    j1j2 <- lookup(mlat,olat)
    i1i2 <- lookup(mlon,olon)
    i1<-i1i2[1]
    i2<-i1i2[2]
    j1<-j1j2[1]
    j2<-j1j2[2]
    I1 <- sub2ind(c(nlon,nlat),i1,j1)
    I2 <- sub2ind(c(nlon,nlat),i2,j1)
    I3 <- sub2ind(c(nlon,nlat),i1,j2)
    I4 <- sub2ind(c(nlon,nlat),i2,j2)
    T1 <- grid.y[I1] - LapseRate*(data.alt[i]-grid.alt[I1])/1000.0
    T2 <- grid.y[I2] - LapseRate*(data.alt[i]-grid.alt[I2])/1000.0
    T3 <- grid.y[I3] - LapseRate*(data.alt[i]-grid.alt[I3])/1000.0
    T4 <- grid.y[I4] - LapseRate*(data.alt[i]-grid.alt[I4])/1000.0
    if (method=='bilinear') {
      y[i] <- intcoef(olon,olat,mlon[i1],mlon[i2],mlat[j1],mlat[j2])%*%as.matrix(c(T1,T2,T3,T4),nrow=4)
    } else if (method =='nearest') {
      h4 <- intcoef.nearest(olon,olat,mlon[i1],mlon[i2],mlat[j1],mlat[j2])
      y[i] <- as.matrix(c(T1,T2,T3,T4),nrow=4)[which(h4==1)[1]]
    } else {
      stop('unknown methods in interpolation')
    }
  }
  
  return(y)
  
}

#' @export
grid2points_test2<-function(grid,data,variable="temperature"){
  grid.grid <- sp::getGridTopology(grid)
  elon<-sp::coordinatevalues(grid.grid)$longitude
  elat<-sp::coordinatevalues(grid.grid)$latitude
  H<-fastgrid::Hmat(elon,elat,coordinates(data))
  y <- as.matrix(H%*%as.matrix(grid@data[,variable]))
  return(y)
}