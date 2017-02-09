!!! Krieging code to speed up computations done in R

!!! R compile:
!!! R CMD SHLIB kriegecode.f90 -lRblas -lRlapack

!!! this version does not use shape information and does not have internal checks

!! kriging predictions as a linear model
!! y = x*beta + eps, eps ~ N(0,b)
!! ypred = xgrid*beta + Bgrid*B^(-1)(y-x*beta)
!! cy, cgrid coordinates of data and the grid
subroutine kriegepred(x,y,b,grid, ypred, cy,cgrid,nobs,npar,ngrid, covpars)
  implicit none
  
  !  include 'lapack_inc.f90'
  integer, parameter :: Double = selected_real_kind(15)
  integer, parameter :: dbl=Double
  integer, parameter :: ncovpars = 2

  integer, intent(in) :: nobs, npar, ngrid
  real(kind=dbl), intent(inout) :: x(nobs,npar), y(nobs,1), b(nobs,nobs), grid(ngrid,npar), ypred(ngrid)
  real(kind=dbl), intent(in) :: cy(nobs,2), cgrid(ngrid,2)
  real(kind=dbl), intent(in) :: covpars(ncovpars)
  
  real(kind=dbl) :: beta(npar,1)
  real(kind=dbl) :: x1(npar), b1(nobs)
    
  real(kind=dbl) :: ypred1(1)
    
  integer :: i

  real(kind=dbl), parameter :: one = 1.0d0, zero = 0.0d0
   
  call chol(b, nobs)               ! replace B with L = chol(B)
  call ldivide(b,x, nobs, npar)    ! replace x with L\x
  call ldivide(b,y, nobs, 1)       ! replace y with L\y
  call lsq(x,y,beta, nobs, npar)   ! beta = x\y
  call resid(x,y,beta, nobs, npar) ! y = y - x*beta
  call ltdivide(b,y, nobs, 1)      ! y = L'\y

  ! calculate prediction over grid, ypred(i) = X1*beta + B1*B\(y-X*beta)
  do i=1, ngrid
     x1 = grid(i,:)
     call corrdist(cgrid(i,:),cy,b1, nobs, covpars)
     ypred1 = matmul(x1,beta) + matmul(b1,y)
     ypred(i) = ypred1(1)
  end do

contains
!!! Utility functions to call LAPACK and BLAS library routines for some matrix computations.
!!! The utilites do not use shape information, so they can be called from R.
!!! Calls: dpotrf, dtrsm, dgels, dgemv
  
  !! correlation calculations
  subroutine corrdist(x1,x,b1, n, covpars)
    implicit none
    real(kind=dbl), intent(in) :: x1(:), x(:,:)
    integer, intent(in) :: n
    real(kind=dbl), intent(out) :: b1(n)
    real(kind=dbl), intent(in) :: covpars(ncovpars)
    
    integer :: i
    real(kind=dbl) :: d
    real(kind=dbl) :: sig2, clen
    
    sig2 = covpars(1)
    clen = covpars(2)
    
    do i=1,n
       d = sqrt((x1(1)-x(i,1))**2 + (x1(2)-x(i,2))**2)
       b1(i) = exp(-d/clen)*sig2
    end do  
  end subroutine corrdist

!! land sea mask
  subroutine corrdist2(x1,x,b1, lsm, n, covpars)
    implicit none
    real(kind=dbl), intent(in) :: x1(:), x(:,:), lsm(:)
    integer, intent(in) :: n
    real(kind=dbl), intent(out) :: b1(n)
    real(kind=dbl), intent(in) :: covpars(ncovpars)
    
    integer :: i
    real(kind=dbl) :: d
    real(kind=dbl) :: sig2, clen
    
    sig2 = covpars(1)
    clen = covpars(2)
    
    do i=1,n
       d = sqrt((x1(1)-x(i,1))**2 + (x1(2)-x(i,2))**2)
       if (lsm(i) .eq. 0.0d0) then
          if (d .eq. 0.0d0) then
             b1(i) = sig2
          else
             b1(i) = 0.0d0
          end if
       else
          b1(i) = exp(-d/clen)*sig2
       end if
    end do  
  end subroutine corrdist2

  
  !! X = chol(X,'lower'), X(1,1) == -1 on error
  subroutine chol(x,n)
    implicit none
    integer, intent(in) :: n
    real(kind=dbl),intent(inout) :: x(n,n)

    integer :: info

    if (n<1) return
    call dpotrf('l',n,x,n,info)

    if (info .ne. 0) then
       x(1,1) = -1.0d0      
       if (n > 1) x(1,2) = info
    end if    
  end subroutine chol

  !! X = L\X, L is lower chol, replaces X
  subroutine ldivide(l,x, n, m)
    implicit none
    integer, intent(in) :: n, m
    real(kind=dbl),intent(in) :: l(n,n)
    real(kind=dbl),intent(inout) :: x(n,m)

    call dtrsm('left','l','n','n',n,m,one,l,n,x,n)
  end subroutine ldivide

  !! X = L'\X, L is lower chol, replaces X
  subroutine ltdivide(l,x, n,m)
    implicit none
    integer, intent(in) :: n, m
    real(kind=dbl),intent(in) :: l(n,n)
    real(kind=dbl),intent(inout) :: x(n,m)

    call dtrsm('left','l','t','n',n,m,one,l,n,x,n)
  end subroutine ltdivide

  !! beta = x\y, general least squares solution, now assumes size(y,2) == 1!
  subroutine lsq(x,y,beta, m, n) ! oops n and m in wrong way!!
    implicit none
    integer, intent(in) :: n, m
    real(kind=dbl), intent(in) :: x(m,n)
    real(kind=dbl), intent(in) :: y(m,1)
    real(kind=dbl), intent(out) :: beta(n)

    integer :: nrhs, lwork, info
    real(kind=dbl), allocatable :: work(:), x2(:,:), y2(:,:)
    
    nrhs = 1 ! size(y,2)
    lwork = min(m,n)+max(1,m,n,nrhs)
    allocate(work(lwork), x2(m,n), y2(m,nrhs))
    y2 = y
    x2 = x
    call dgels('n', m, n, nrhs, x2, m, y2, m, work, lwork, info)
    if (info .ne. 0) then
!       write(*,*) 'dgels failed'
!       write(*,*) 'info = ', info
    endif
    beta = y2(1:n,1)
    deallocate(work,x2,y2)
  end subroutine lsq

  !! residuals y = y-x*beta, replaces y
  subroutine resid(x,y,beta, n,m)
    implicit none
    integer, intent(in) :: m,n
    real(kind=dbl), intent(in) :: x(n,m)
    real(kind=dbl), intent(inout) :: y(n,1)
    real(kind=dbl), intent(in) :: beta(m)
    call dgemv('n',n,m,-one,x,n,beta,1,one,y,1)
  end subroutine resid
    
end subroutine kriegepred

subroutine obsoper(mlat,mlon,obs,nlat,nlon,nobs,iind,jind,h)
  implicit none

  integer, parameter :: Double = selected_real_kind(15)
  integer, parameter :: dbl=Double

  integer, intent(in) :: nlat, nlon, nobs
  real(kind=dbl), intent(in) :: mlat(nlat), mlon(nlon), obs(nobs,2)

  integer, intent(out) :: iind(4*nobs), jind(4*nobs)
  real(kind=dbl), intent(out) :: h(4*nobs)

  integer :: nmod, i, ii, j, i1i2(2),j1j2(2), i1,i2,j1,j2, ii4(4)
  real(kind=dbl) :: olat, olon, h4(4)

  nmod = nlat*nlon
!  # sparse matrix
!  H<-Matrix(0,nrow=nobs,ncol=nmod,sparse=TRUE)

  j = 1
  do i = 1,nobs
    olat = obs(i,1)
    olon = obs(i,2)
    i1i2 = lookup(mlat,olat)
    j1j2 = lookup(mlon,olon)
    i1 = i1i2(1)
    i2 = i1i2(2)
    j1 = j1j2(1)
    j2 = j1j2(2)
    II4(1) = sub2ind(nlat,i1,j1)
    II4(2) = sub2ind(nlat,i2,j1)
    II4(3) = sub2ind(nlat,i1,j2)
    II4(4) = sub2ind(nlat,i2,j2)
    h4 = intcoef(olat,olon,mlat(i1),mlat(i2),mlon(j1),mlon(j2))
    do ii=1,4
       h(j) = h4(ii)
       iind(j) = i
       jind(j) = II4(ii)
       j = j + 1
    end do
 end do
 
contains
  
  function intcoef(x1,x2,x1l,x1u,x2l,x2u) result(h)
    implicit none
    real(kind=dbl), intent(in) :: x1,x2,x1l,x1u,x2l,x2u
    real(kind=dbl) :: h(4)
    
    real(kind=dbl) :: X(4,4), x11, x22

    X = reshape((/1.0,-1.0,-1.0,1.0,0.0,1.0,0.0,-1.0,0.0,0.0,1.0,-1.0,0.0,0.0,0.0,1.0/),(/4,4/))
    x11 = (x1-x1l)/(x1u-x1l)
    x22 = (x2-x2l)/(x2u-x2l)
    h = reshape(matmul(reshape((/1.0d0,x11,x22,x11*x22/),(/1,4/)),X),(/4/))
  end function intcoef

  function sub2ind(nrc,irow,icol) result(i)
    implicit none
    integer, intent(in) :: nrc, irow, icol
    integer :: i
    i = (icol-1)*nrc + irow
  end function sub2ind

  function lookup(x,xin) result(ii)
    implicit none
    real(kind=dbl) :: x(:), xin
    integer :: ii(2)
    integer :: n, il, iu, im

    n = size(x,1)
    !!  find location by bisection
    if (x(2)-x(1)>0.0) then
       il=1
       iu=n
       do while (iu-il>1)
          im = floor((iu+il)/2.0)
          if (xin>x(im)) then
             il=im
          else
             iu=im
          end if
       end do
       ii = (/il,iu/)
    else
       il=n
       iu=1
       do while (il-iu>1)
          im = floor((iu+il)/2.0)
          if (xin>x(im)) then
             il=im
          else
             iu=im
          end if
       end do
       ii = (/iu,il/)
    end if
  end function lookup

end subroutine obsoper
