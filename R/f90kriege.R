### Fast Krieging by calling fortran code from R
### INTENDED FOR MOS temperature data in predefined format
### but might work for other purposes, too.
###
### compilation of the shared code by a command like:
### R CMD SHLIB kriegecode.f90 -lRblas -lRlapack
###
### or preferably install as package:
### R CMD INSTALL fastgrid_n.n.tar.gz

### marko.laine@fmi.fi

## TODO: error processing!!!

#' fastgrid: faster Kriging by some f90 code
#'
#' These are some functions here to help with gridding point data to a regular grid.
#' 
#' @section functions:
#' See \code{\link{fastkriege}}.
#'
#'#' @seealso
#' \code{\link[gstat]{gstat}} \code{\link[MOSfieldutils]{MOSfieldutils}}
#' 
#' @docType package
#' @name fastgrid
NULL

# methods package needed if run by Rscript
# but seems to need library(methods) in the file.R too (at least for the older version of Matrix)
# import sp for spatial data and Matrix for sparse matrix used in observation opeerator
#' @import sp Matrix
NULL

## This loads the code to R memory
# if (!is.loaded("kriegecode")) {
#    dyn.load("kriegecode.so")
#}


#' @export
fastkriege_old <- function(trend_model = temperature ~ -1, data, grid, cov.pars, 
                       lsm=NULL,lsmy=NULL, alt=NULL, alty=NULL, altlen=200.0,
                       bg=NULL, variable="temperature" ) {
  ## build input matrices

  trend_model_noy <- trend_model[-2]
  grid.variables <- all.vars(trend_model_noy)
  
  if (!is.null(bg)) {
    ## map from from bg grid to station locations
    ## need the grid for obs operator, FIX me
    s<-sp::summary(grid)$grid
    elon<-seq.int(from=s[1,1],by=s[1,2],len=s[1,3])
    elat<-seq.int(from=s[2,1]+s[2,2]*(s[2,3]-1) ,by=-s[2,2],len=s[2,3])
    H<-f90Hmat(elon,elat,coordinates(data))
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
  
  ## Kriging by fortran code (now 3 separate code for different situations)
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
    B <- B*exp(-as.matrix((dist(alty)/altlen)))
    ypred<-f90kriege3(X,y,B,predgrid,lsm,lsmy,alt,alty,cy,cgrid,c(cov.pars[1],cov.pars[2],altlen))
  }
  t2<-proc.time()-t1
  
  ypredgrid<-double(length=length(igrid))*NA
  ypredgrid[igrid]<-ypred
  
  ## add bg field to the results
  if (is.null(bg))
    ypred2<-data.frame(temperature=ypredgrid)
  else
    ypred2<-data.frame(temperature=ypredgrid+bg@data[,variable])
  ## change the name according to the ´variable´
  names(ypred2) <- c(variable)
  coordinates(ypred2)<-coordinates(grid)
  proj4string(ypred2)<-proj4string(grid)
  gridded(ypred2)<-TRUE
  fullgrid(ypred2) <- TRUE
  
  if (!is.null(bg))
    ypred2$diff <- as.vector(ypred2@data[,variable] - bg@data[,variable])
  
  attr(ypred2,"elapsed") <- t2
  attr(ypred2,'failed') <- attr(ypred,'failed')
  
  return(ypred2)
}

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
                       bg=NULL, variable="temperature", LapseRate = 0.0,
                       method='bilinear') {
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

    if (LapseRate == 0.0) {
      mu <- grid2points(bg, data, variable=variable, method=method)
    } else {
      mu <- grid2points_lapserate(bg,data,modelgrid = grid,
                                  LapseRate = LapseRate, variable=variable,
                                  method=method)
    }
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
  
#  ypredgrid<-double(length=nrow(grid))*NA
  ypredgrid<-double(length=length(igrid))*NA
  ypredgrid[igrid]<-ypred
  
  ## add bg field to the results
  if (is.null(bg))
    ypred2<-data.frame(temperature=ypredgrid)
  else
    ypred2<-data.frame(temperature=ypredgrid+bg@data[,variable])
  ## change the name according to the ´variable´
#   names(ypred2) <- variable
  names(ypred2) <- c(variable)
  coordinates(ypred2)<-coordinates(grid)
  proj4string(ypred2)<-proj4string(grid)
  gridded(ypred2)<-TRUE
  fullgrid(ypred2) <- TRUE # ok? makes the variable smaller
  
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

