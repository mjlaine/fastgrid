# Some grid utilities

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
### OBS changed lat,lon to lon,lat in 2017-12-1!!!
#' Linear observation operator
#' 
#' @param mlon,mlat model longitude (x) and latitude (y) coordinates
#' @param obs matrix with two columns of data x and y coordinates
#' @param method interpolation method, default 'bilinear'
#' 
#' @return 
#' 
#' Sparse Matrix on size nrows(obs) times length(mlon)*length(mlat),
#' each row has at most 4 non-zero elements
#' 
#' @export
Hmat <- function(mlat,mlon,obs,method='bilinear') {
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

    if (method=='bilinear') {
      H[i,c(I1,I2,I3,I4)] <- intcoef(olat,olon,mlat[i1],mlat[i2],mlon[j1],mlon[j2])
    } else if (method =='nearest') {
      H[i,c(I1,I2,I3,I4)] <- intcoef.nearest(olat,olon,mlat[i1],mlat[i2],mlon[j1],mlon[j2])
    } else {
      stop('unknown methods in Hmat')
    }
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
  # order is ll,ul,lu,uu
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




# utility for SpatialPixelsDataFrame coordinates
#' Grid coordinates
#'
#' Get longitude or latitude coordinate values of a regular spatial grid.
#' Uses \code{\link[sp]{getGridTopology}} and  \code{\link[sp]{coordinatevalues}} from the \code{\link[sp]{sp}} package.
#'
#' @param x spatial grid
#' @return distinct coordinate values
#'
#' @examples
#'
#' gridlat(MOSgriddata())
#'
#' @export
gridlon <- function(x) {
  sp::coordinatevalues(sp::getGridTopology(x))$longitude
}
#' @rdname gridlon
#' @export
gridlat <- function(x) {
  sp::coordinatevalues(sp::getGridTopology(x))$latitude
}

# this does not use coordinate names
#    s<-sp::summary(x)$grid
#    elon<-seq.int(from=s[1,1],by=s[1,2],len=s[1,3])
#    elat<-seq.int(from=s[2,1]+s[2,2]*(s[2,3]-1) ,by=-s[2,2],len=s[2,3])

## map grid values to points (uses H from fastgrid)
#' Map grid to poinits
#' 
#' @param grid SpatialGridDataFrame
#' @param data SpatialPointsDataFrame
#' @param variable variable to be mapped, default \code{temperature}
#' 
#' @export
grid2points<-function(grid,data,variable="temperature",method='bilinear'){
  grid.grid <- sp::getGridTopology(grid)
  elon<-sp::coordinatevalues(grid.grid)$longitude
  elat<-sp::coordinatevalues(grid.grid)$latitude
  
  if (method=='bilinear') {
    H<-f90Hmat(elon,elat,coordinates(data))
    y <- as.matrix(H%*%as.matrix(grid@data[,variable]))
  } else if (method =='nearest') {
    H<-Hmat(elon,elat,coordinates(data),method=method)
    y <- as.matrix(H%*%as.matrix(grid@data[,variable]))
  } else {
    stop('unknown methods in interpolation')
  }
  
  return(y)
}

## map grid values to points (there is one in MOSfieldutils, too)
#' @rdname grid2points
#' @export
grid2points_lapserate <- function(grid,data,variable="temperature", 
                                  LapseRate=MOSget('LapseRate'),
                                  modelgrid=NULL, method='bilinear') {
  grid.grid <- sp::getGridTopology(grid)
  # coordinate names fixed here!!
  mlon<-sp::coordinatevalues(grid.grid)$longitude
  mlat<-sp::coordinatevalues(grid.grid)$latitude
  
  nlat <- length(mlat)
  nlon <- length(mlon)
  nmod <- nlat*nlon
  
  # needs MOSfieldutils here  
  if (is.null(modelgrid)) {
    # data("KriegeData", package = MOS.options$pkg, envir = parent.frame())
    modelgrid <- MOSfieldutils::MOSgriddata()
  }
  
  grid.y <- as.matrix(grid@data[,variable])
  data.coord <- sp::coordinates(data)
  if (LapseRate != 0.0) {
    #  grid.alt <- grid$elevation
    grid.alt <- modelgrid$elevation
    data.alt <- data$elevation
  }
  
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
    if (LapseRate != 0.0) {
      T1 <- grid.y[I1] - LapseRate*(data.alt[i]-grid.alt[I1])/1000.0
      T2 <- grid.y[I2] - LapseRate*(data.alt[i]-grid.alt[I2])/1000.0
      T3 <- grid.y[I3] - LapseRate*(data.alt[i]-grid.alt[I3])/1000.0
      T4 <- grid.y[I4] - LapseRate*(data.alt[i]-grid.alt[I4])/1000.0
    } else {
      T1 <- grid.y[I1]
      T2 <- grid.y[I2]
      T3 <- grid.y[I3]
      T4 <- grid.y[I4]
    }
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


#' Generate uniform rectangular longitude-latitude-based SpatialGrid
#' 
#' @param lon,lat longitude and latitude grid center points
#' @param lonmin,lonmax,londx longitude min, max and cell size, used if \code{lon} not given.
#' @param latmin,latmax,latdx latitude min, max and cell size, used if \code{lat} not given.
#' 
#' @examples  
#' 
#' sp::getGridTopology(makelonlatgrid(10,20,-10,10))
#' 
#' @export
makelonlatgrid <- function(lonmin, lonmax, latmin, latmax, londx=1.0, latdx=1.0,
                     lon=NULL, lat=NULL) {
  
  if (!is.null(lon)) {
    lonmin <- min(lon)
    lonmax <- max(lon)
    londx <- abs(lon(2)-lon(1))
  }
  if (!is.null(lat)) {
    latmin <- min(lat)
    latmax <- max(lat)
    latdx <- abs(lat(2)-lat(1))
  }
  modelgrid <- sp::SpatialGrid(sp::GridTopology(cellcentre.offset=c(lonmin,latmin),
                                                cellsize=c(londx,latdx),
                                                cells.dim=c(ceiling((lonmax-lonmin)/londx)+1,
                                                            ceiling((latmax-latmin)/latdx)+1)),
                               proj4string=sp::CRS("+init=epsg:4326"))
  sp::coordnames(modelgrid)<-c('longitude','latitude')
  return(modelgrid)
  
}


# utility to points2grid
points2grid0 <- function(pointdata,modelgrid,variable="y",
                         trend_model=NULL,bg=NULL,cov.pars) {
  if (is.null(trend_model)) trend_model <- as.formula(paste(variable,'~1'))
  pointdata <- pointdata[complete.cases(pointdata@data[,variable]),]
  if (length(pointdata)<4)  stop('too few points')
  Xgrid <- fastgrid::fastkriege(trend_model=trend_model, data=pointdata, grid=modelgrid, cov.pars=cov.pars,
                                bg=bg,lsm=NULL,lsmy=NULL,alt=NULL,alty=NULL,variable=variable)
  return(Xgrid)
}

# can not use points2grid as it is already used by sp:points2grid

#' Grid point data to spatial grid
#' 
#' This is an easier to use version of \link{fastkriege}
#' 
#' @param y data matrix of point data to be gridded, alternative to giving \code{variable}
#' @param pointdata SpatialPointsDataFrame of the spatial points
#' @param modelgrid SpatialGrid or SpatialGridDataFrame of the output grid
#' @param variable name of the variable in \code{pointdata} to be gridded
#' @param cov.par Kriging parameters, sigma^2, phi, tau^2.
#' 
#' @examples  
#' 
#' lonmin <- 19.0; lonmax <- 33.0
#' latmin <- 59.0; latmax <- 71.5
#' newgrid <- makelonlatgrid(lonmin=lonmin,lonmax=lonmax,latmin=latmin,latmax=latmax,londx=0.1,latdx=0.1)
#' nobs <- 50
#' obs <- data.frame(longitude=runif(nobs,lonmin,lonmax),
#'                  latitude=runif(nobs,latmin,latmax),
#'                  temperature=rnorm(nobs,mean=10,sd=4))
#' sp::coordinates(obs) <- c('longitude','latitude')
#' out <- pointgridding(pointdata=obs, modelgrid = newgrid, variable="temperature",cov.pars = c(0.5^2,2.0,0.0))
#' MOSplotting::MOS_plot_field(out,stations=obs,cmin=0,cmax=20)
#' 
#' @export
pointgridding <- function(y=NULL,pointdata,modelgrid,variable="y",trend_model=NULL,
                          priorfield=FALSE,cov.pars) {
  
  if (priorfield) {
    bg <- modelgrid
  } else {
    bg <- NULL
  }
  if (is.null(y)) {
    return(points2grid0(pointdata=pointdata, modelgrid = modelgrid, bg=bg, variable=variable, trend_model, cov.pars = cov.pars))
  }
  
  if (is.null(trend_model)) trend_model <- as.formula(paste(variable,'~1'))
  
  # else use each column in matrix y
  y <- as.matrix(y)
  pointdata$y <- y[,1]
  out <- points2grid0(pointdata = pointdata, modelgrid = modelgrid, bg=bg, variable="y",trend_model, cov.pars = cov.pars)
  
  l <- dim(y)[2]
  if (l>1) {
    for (i in 2:l) {
      pointdata$y <- y[,i]
      out.1 <- points2grid0(pointdata = pointdata, modelgrid = modelgrid, bg=bg, variable="y", trend_model, cov.pars = cov.pars)
      ii <- dim(out@data)[2]+1
      out@data[,ii] <- out.1@data[,"y"]
    }
  }
  return(out)
}

