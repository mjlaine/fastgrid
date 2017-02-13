### Fast Krieging by calling fortran Krieging code from R
### This only works with MOS temperature data in predefined format
###
### compilation of the shared code by a command like:
### R CMD SHLIB kriegecode.f90 -lRblas -lRlapack
###
### or install as package:
### R CMD INSTALL fastgrid_n.n.tar.gz

### marko.laine@fmi.fi

## This loads the code to R memory
# if (!is.loaded("kriegecode")) {
#    dyn.load("kriegecode.so")
#}

#require(Matrix) # for sparseMatrix
#require(sp) # for spatial objects

## cover function for fortran kriging code in kriegecode.f90
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
f90kriege2 <- function(x,y,b,grid,lsm,cy,cgrid,covpars) {
  nobs <- nrow(x)
  npar <- ncol(x)
  ngrid <- nrow(grid)
  
  ## There is very little checks in the Fortran code, so check the sanity of the inputs here 
  if (nrow(y) != nobs | ncol(y) != 1 | nrow(cy) != nobs
      | ncol(grid) != npar | nrow(cgrid) != ngrid | length(lsm) != ngrid
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



### cover function for fortran code for observation operator
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

  return(sparseMatrix(i=f90h$iind,j=f90h$jind,x=f90h$h,dims=c(nobs,nmod)))
}

## to replace gstat kriege command for prediction over a grid
fastkriege <- function(trend_model, data, grid, cov.pars, lsm=NULL, bg=NULL,variable="temperature" ) {
    ## build input matrices
    ## assumes data and bg have "longitude", "latitude", "temperature", and data and grid also trend model variables

    trend_model_noy <- trend_model[-2]
    grid.variables <- all.vars(trend_model_noy)
  
    if (!is.null(bg)) {
      ## map from from bg grid to station locations
      ## need the grid for obs operator, FIX me
      s<-summary(grid)$grid
      elon<-seq.int(from=s[1,1],by=s[1,2],len=s[1,3])
      elat<-seq.int(from=s[2,1]+s[2,2]*(s[2,3]-1) ,by=-s[2,2],len=s[2,3])
      H<-f90Hmat(elon,elat,cbind(data$longitude,data$latitude))
#      mu <- as.matrix(H%*%bg$temperature)
      mu <- as.matrix(H%*%as.matrix(bg@data[variable]))
    }
    else {
        mu <- 0.0
    }

    ## non missing of all trend model variables
    igrid <- complete.cases(as.data.frame(grid)[,grid.variables])
  
    B <- buildcovmat(coords=coordinates(data), cov.model = "exp", cov.pars=cov.pars)
    X <- model.matrix(trend_model_noy, data=data)
#    y <- as.matrix(data$temperature-mu)
    y <- as.matrix(data@data[variable]-mu)
    predgrid<-model.matrix(trend_model_noy,data=grid)
    cy <- as.matrix(coordinates(data))
    cgrid <- as.matrix(coordinates(grid))[igrid,]
    
    ## Kriging by fortran code
    t1<-proc.time()
    if (is.null(lsm))
      ypred<-f90kriege(X,y,B,predgrid,cy,cgrid,cov.pars)
    else
      ypred<-f90kriege2(X,y,B,predgrid,lsm,cy,cgrid,cov.pars)
    end
    t2<-proc.time()-t1

    ypredgrid<-double(length=nrow(grid))*NA
    ypredgrid[igrid]<-ypred
  
    ## add bg field to the results
    if (is.null(bg))
      ypred2<-data.frame(temperature=ypredgrid)
    else
  #      ypred2<-data.frame(temperature=ypredgrid+bg$temperature)
      ypred2<-data.frame(temperature=ypredgrid+bg@data[variable])
    names(ypred2) <- c(variable)
    coordinates(ypred2)<-coordinates(grid)
    proj4string(ypred2)<-proj4string(grid)
    gridded(ypred2)<-TRUE
    return(ypred2)
}


### Observation operator using bi-linear interpolation and bisection
### Returns sparse matrix
### See also: f90Hmat, which uses fortran code
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
                    cov.pars = stop("no cov.pars argument"), dist.method = "euclidian", ...) {
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


