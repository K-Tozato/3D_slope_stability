!----------------------------------
module coord 
  integer :: nx, ny
  double precision :: dx,dy
  double precision, allocatable :: xx(:), yy(:), gl(:,:)
  double precision, allocatable :: xxn(:), yyn(:), gln(:,:)
end module

module inp_data
  double precision, allocatable :: gd(:,:), gs(:,:), pd(:,:) 
  double precision, allocatable :: ps(:,:), cd(:,:), cs(:,:)
  double precision :: sra(2)
end module

module ssa
  integer, allocatable :: nsl(:,:), nep(:,:), nell(:,:)
  double precision, allocatable :: uu(:,:), vv(:,:), ww(:,:)
  double precision, allocatable :: pell(:,:)
  double precision, allocatable :: fsp(:), zvp(:)
end module

module solve_data
  double precision, allocatable :: fs(:,:), zz(:,:), zvol(:,:)
end module

!==================================
program hovland
use solve_data, only : fs
implicit none
integer :: nn, ne, i, nf, val(8)
integer, allocatable :: ngw(:)
character(10) :: ct(3)

open(11,file='./input/gwater_case.txt',status='old')
read(11,*) ne
read(11,*) nn
close(11)
nf = nn - ne + 1
allocate(ngw(nf))
do i = 1,nf
  ngw(i) = ne + i - 1
end do

call datain
call angle(nn)
call ell_size(ne)

do i = 1, nf
  write(*,*) '==============='
  write(*,'(a,i4,a)') ' | Case : ', ngw(i), ' |----------'
  write(*,*) '==============='
  call date_and_time(ct(1),ct(2),ct(3),val)
  write(*,'(a,i4,a,i2,a,i2, i6,a,i2,a,i2)') 'Start : ',&
  val(1),'/',val(2),'/',val(3), val(5),':',val(6),':',val(7)
  
  call hovland_calc(nn,ne,ngw(i))

  write(*,'(a,f8.4)') 'Minimum Fs : ', minval(fs)
  call date_and_time(ct(1),ct(2),ct(3),val)
  write(*,'(a,i4,a,i2,a,i2, i6,a,i2,a,i2)') ' End  : ', &
  val(1),'/',val(2),'/',val(3), val(5),':',val(6),':',val(7)
  write(*,*) '-------------------------'
end do
end program
!==================================

! - - - - - - - - - - - - - 
subroutine datain
use coord
use ssa
use inp_data 
use solve_data
implicit none
integer :: i, j, k, nn
double precision :: dum(2)
double precision, parameter :: pi = 4.0d0*atan(1.d0)

open(11,file='./input/num_node.txt',status='old')
read(11,*) nx
read(11,*) ny
close(11)
allocate(xx(nx-1),yy(ny-1),gl(nx-1,ny-1))
allocate(xxn(nx),yyn(ny),gln(nx,ny))
allocate(fs(nx-1,ny-1),zz(nx-1,ny-1),zvol(nx-1,ny-1))
allocate(cd(nx-1,ny-1), pd(nx-1,ny-1), gd(nx-1,ny-1)) 
allocate(cs(nx-1,ny-1), ps(nx-1,ny-1), gs(nx-1,ny-1))

fs(:,:) = 10.d0 ; zvol(:,:) = 0.d0

open(11,file='./input/coordinate.txt',status='old')
nn = nx*ny
i = 0  ; j = 1
do k = 1,nn
  i = i + 1
  read(11,*) xxn(i), yyn(j), gln(i,j)
  if(i .eq. nx) then
    i = 0 ; j = j + 1
  end if
end do
close(11)

dx = xxn(2) - xxn(1)
dy = yyn(2) - yyn(1)

open(11,file='./input/parameter_slope.txt',status='old')
do j = 1,ny-1
  do i = 1, nx-1
    xx(i)   = 0.50d0 * (xxn(i+1) + xxn(i))
    yy(j)   = 0.50d0 * (yyn(j+1) + yyn(j))
    gl(i,j) = 0.25d0 * (gln(i,j) + gln(i+1,j) + gln(i,j+1) + gln(i+1,j+1))  
    read(11,*) cd(i,j), pd(i,j), gd(i,j), cs(i,j), ps(i,j), gs(i,j)
  end do
end do
close(11)
ps(:,:) = ps(:,:)/180.d0*pi ; pd(:,:) = pd(:,:)/180.d0*pi

open(11,file='./input/infsupdip.dat',status='old')
read(11,*) sra(1), sra(2)
close(11)
sra(1) = sra(1)/180.d0*pi ; sra(2) = sra(2)/180.d0*pi

end subroutine 
! - - - - - - - - - - - - - 
subroutine angle(nn)
use coord
use ssa
use inp_data, only : sra
implicit none
integer :: i, j, simple, is, js
integer, intent(out) :: nn
double precision :: dzx, dzy, slp, tth, dd

!write(*,*) "Simple slope ?"
!read(*,*) simple
simple = 0
if(simple==1) then
  open(11,file='./input/simple_slope.txt',status='old')
  read(11,*) is
  read(11,*) js
  close(11)
end if

nn = 0
do j = 1, ny-1
  do i = 1, nx-1 
    dzx = (gln(i+1,j+1)+gln(i+1,j)-gln(i,j+1)-gln(i,j))*0.5d0/dx
    dzy = (gln(i+1,j+1)+gln(i,j+1)-gln(i+1,j)-gln(i,j))*0.5d0/dy
    tth = dsqrt(dzx**2.d0 + dzy**2.d0)
    slp = datan(tth)
    if(slp .ge. sra(1) .and. slp .le. sra(2)) then
      nn = nn + 1 ; cycle
    end if
  end do
end do
allocate(nsl(nn,2),uu(nn,3),vv(nn,3),ww(nn,3))

nn = 0
do j = 1, ny-1
  do i = 1, nx-1 
    if(simple == 1) then
      if(i .ne. is .or. j .ne. js) cycle
    end if

    dzx = (gln(i+1,j+1)+gln(i+1,j)-gln(i,j+1)-gln(i,j))*0.5d0/dx
    dzy = (gln(i+1,j+1)+gln(i,j+1)-gln(i+1,j)-gln(i,j))*0.5d0/dy
    tth = dsqrt(dzx**2.d0 + dzy**2.d0)
    slp = datan(tth)

    if(slp .ge. sra(1) .and. slp .le. sra(2)) then
      nn = nn + 1
      nsl(nn,1) = i ; nsl(nn,2) = j
      
      !dzx = (gln(i+1,j+1)+gln(i+1,j)-gln(i,j+1)-gln(i,j))*0.5d0
      !dzy = (gln(i+1,j+1)+gln(i,j+1)-gln(i+1,j)-gln(i,j))*0.5d0
      dd  = dsqrt(dzx**2.d0 + dzy**2.d0 +(dzx**2.d0+dzy**2.d0)**2.d0)
      uu(nn,1) = -  dzx / dd
      uu(nn,2) = -  dzy / dd
      uu(nn,3) = - (dzx**2.d0+dzy**2.d0) / dd
      
      ww(nn,1) =   dzx / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
      ww(nn,2) =   dzy / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
      ww(nn,3) = -1.d0 / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
      
      vv(nn,1) = uu(nn,2)*ww(nn,3) - uu(nn,3)*ww(nn,2)
      vv(nn,2) = uu(nn,3)*ww(nn,1) - uu(nn,1)*ww(nn,3)
      vv(nn,3) = uu(nn,1)*ww(nn,2) - uu(nn,2)*ww(nn,1)
    end if

    if(simple == 1 .and. nn==1) exit
  end do 
  if(simple == 1 .and. nn==1) exit
end do

write(*,*) 'Number of grid: ', (nx-1)*(ny-1)
write(*,*) 'Number of evaluation points:', nn

end subroutine
! - - - - - - - - - - - - - 
subroutine ell_size(ne)
use coord, only : dx
use ssa, only : pell, nep, fsp, zvp, nell
implicit none
integer :: i,j,k,l,ii,jj,kk,ll,n
integer, intent(out) :: ne
double precision, allocatable:: ru(:), rv(:), rw(:), ctr(:)


open(11,file='./input/ru.txt',status='old')
open(12,file='./input/rv.txt',status='old')
open(13,file='./input/rw.txt',status='old')
open(14,file='./input/center.txt',status='old')
read(11,*) ii
read(12,*) jj
read(13,*) kk
read(14,*) ll
allocate(ru(ii),rv(jj),rw(kk),ctr(ll))
do i = 1,ii
  read(11,*) ru(i)
end do
close(11)
do i = 1,jj
  read(12,*) rv(i)
end do
close(12)

do i = 1,kk
  read(13,*) rw(i)
end do
close(13)

do i = 1,ll
  read(14,*) ctr(i)
end do
close(14)

allocate(pell(ii*jj*kk*ll,4))
n = 0
do i = 1,ii
  do j = 1,jj
    do k = 1,kk
      do l = 1,ll
        n = n + 1
        pell(n,1) = ru(i)
        pell(n,2) = rv(j)*pell(n,1)
        pell(n,3) = rw(k)*pell(n,1)
        pell(n,4) = ctr(l)
      end do
    end do
  end do
end do
write(*,*) 'Number of ellipsoid shape:', n
ne = n

n = int((maxval(pell(:,1))/dx*2.d0)**2.d0*0.8d0)
allocate(nep(n,2),fsp(n),zvp(n),nell(n,2))
end subroutine
! - - - - - - - - - - - - - 
subroutine hovland_calc(nn,ne,nf)
use coord
use ssa
use inp_data, only: gd, gs, pd, ps, cd, cs
use solve_data
implicit none
integer, intent(in) :: nn, ne, nf
integer :: i, j, k, kk, l, ll, i1, j1, n0, posi(4,2), mm(3), n, n1
double precision :: tmat(3,3), tmatt(3,3), mat(3,3)
double precision :: ee(3), ctr(3), xyz(3), val, tij(3), nij(3), gv(3)
double precision :: xp(4), yp(4), zp(4), rb(3), rg(3)
double precision :: ff(2), a, b, c, dzx, dzy, wij, uij, area, phi
double precision, parameter :: gg = 9.81d0, rw = 1000.d0
character(100) :: fname, fname1

fname1 = './output/water-depth-level_'
write(fname,'(a,i4.4,a4)') trim(adjustl(fname1)), nf, '.txt'
open(11,file=fname, status='old')
i = 0 ; j = 1
do k = 1,(nx-1)*(ny-1)
  i = i + 1
  read(11,*) ee(1), ee(2), zz(i,j), ee(3)
  if(i == nx-1) then 
    i = 0 ; j = j + 1
  end if
end do
close(11)
gv = (/0.d0,0.d0,-1.d0/) ; mm  = (/1,2,3/)
posi(:,:) = 0 ; posi(2,1) = 1 ; posi(3,2) = 1 ; posi(4,:) = 1

!$OMP parallel
!$OMP do private(fsp,zvp,nep,nell,i,j,k,kk,l,ll,i1,j1,n0,tmat,tmatt,mat,ee,ctr,xyz,val,tij,nij,xp,yp,zp,rb,rg) &
!$OMP & private(ff,a,b,c,dzx,dzy,wij,uij,area,phi), reduction(MIN:fs), reduction(MAX:zvol)

do n = 1,nn
  do n1 = 1,ne

!!-------------------------------------------------------------------
i1 = nsl(n,1) ; j1 = nsl(n,2)

ctr(1) = xx(i1)    - pell(n1,4)*pell(n1,3)*ww(n,1)
ctr(2) = yy(j1)    - pell(n1,4)*pell(n1,3)*ww(n,2)
ctr(3) = gl(i1,j1) - pell(n1,4)*pell(n1,3)*ww(n,3)

nep = 0 ; fsp = 10.d0 ; zvp = 0.d0 
ee = (/1.d0,0.d0,0.d0/)
tmat(1,1) = dot_product(ee,uu(n,:))
tmat(2,1) = dot_product(ee,vv(n,:))
tmat(3,1) = dot_product(ee,ww(n,:))
ee = (/0.d0,1.d0,0.d0/)
tmat(1,2) = dot_product(ee,uu(n,:))
tmat(2,2) = dot_product(ee,vv(n,:))
tmat(3,2) = dot_product(ee,ww(n,:))
ee = (/0.d0,0.d0,1.d0/)
tmat(1,3) = dot_product(ee,uu(n,:))
tmat(2,3) = dot_product(ee,vv(n,:))
tmat(3,3) = dot_product(ee,ww(n,:))

tmatt(1,:) = tmat(:,1) ; tmatt(2,:) = tmat(:,2) ; tmatt(3,:) = tmat(:,3)
mat(:,:) = 0.d0
mat(1,1) = 1.d0/(pell(n1,1)**2.d0)
mat(2,2) = 1.d0/(pell(n1,2)**2.d0)
mat(3,3) = 1.d0/(pell(n1,3)**2.d0)
mat = matmul(tmatt,mat)
mat = matmul(mat,tmat)

n0 = 0 ; ll = 0
do j = max(1,j1-int(pell(n1,1)/dy)-2), min(ny-1,j1+int(pell(n1,1)/dy)+2)
  do i = max(1,i1-int(pell(n1,1)/dx)-2), min(nx-1,i1+int(pell(n1,1)/dx)+2)
    l = 0
    do k = 1,4 
      xyz(1) = xxn(i+posi(k,1)) - ctr(1)
      xyz(2) = yyn(j+posi(k,2)) - ctr(2)
      xyz(3) = gln(i+posi(k,1), j+posi(k,2)) - ctr(3) 
      ee(1)  = xyz(1)*mat(1,1)+xyz(2)*mat(2,1)+xyz(3)*mat(3,1)
      ee(2)  = xyz(1)*mat(1,2)+xyz(2)*mat(2,2)+xyz(3)*mat(3,2)
      ee(3)  = xyz(1)*mat(1,3)+xyz(2)*mat(2,3)+xyz(3)*mat(3,3)
      val    = dot_product(ee,xyz)

      if(val<=1.d0) l = l + 1 
      if(val> 1.d0) kk = k
    end do
  
    if(l==3) then
      n0 = n0 + 1 ; nep(n0,1) = i ; nep(n0,2) = j
      nell(n0,1) = 3 ; nell(n0,2) = kk
    end if
    if(l==4) then
      n0 = n0 + 1 ; nep(n0,1) = i ; nep(n0,2) = j
      nell(n0,1) = 4 ; nell(n0,2) = 0
    end if
    if(i==1 .or. i==nx-1 .or. j==1 .or. j==ny-1) then
      if(l==4) ll = ll + 1
    end if
  end do
end do

ff = 0.d0 ; zvp = 0.d0
do k = 1,n0
  if(ll>1 .or. n0<10) exit
  i = nep(k,1) ; j = nep(k,2)
  xp = 0.d0 ; yp = 0.d0 ; zp = 0.d0
  do l = 1, 4
    xyz(1) = xxn(i+posi(l,1)) - ctr(1)
    xyz(2) = yyn(j+posi(l,2)) - ctr(2)
    xyz(3) = gln(i+posi(l,1), j+posi(l,2))
    a = mat(3,3) 
    b = (mat(3,1)+mat(1,3))*xyz(1) + (mat(2,3)+mat(3,2))*xyz(2)
    c = mat(1,1)*(xyz(1)**2.d0) + (mat(2,1)+mat(1,2))*xyz(1)*xyz(2) + mat(2,2)*(xyz(2)**2.d0) -1.d0
    zp(l) = ctr(3) + (- b - dsqrt(b**2.d0 - 4.d0 * a * c)) / (2.d0 * a)
    zp(l) = xyz(3) - zp(l)
    if(l==nell(k,2)) zp(l) = 0.d0
    xp(l) = xyz(1) + ctr(1) ; yp(l) = xyz(2) + ctr(2)
    zvp(k) = zvp(k) + zp(l)/dble(nell(k,1))
  end do

  val = zvp(k) * dx * dy

  if(zz(i,j) .ge. zvp(k)) then
    wij = gs(i,j) * val
    uij = zvp(k) * gg * rw * 0.001d0
    c = cs(i,j) ; phi = ps(i,j)
  else
    a = zz(i,j)/zvp(k)
    wij = (gs(i,j)*a + gd(i,j)*(1.d0-a)) * val
    uij = 0.d0
    c = cd(i,j) ; phi = pd(i,j)
  end if
 
  do l = 1,4
    zp(l) = - zp(l) + gln(i+posi(l,1),j+posi(l,2)) 
  end do

  dzx = (zp(4) + zp(2) - zp(3) - zp(1))*0.5d0/dx
  dzy = (zp(4) + zp(3) - zp(2) - zp(1))*0.5d0/dy

  nij(1) =    dzx / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
  nij(2) =    dzy / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
  nij(3) = - 1.d0 / dsqrt(dzx**2.d0 + dzy**2.d0 + 1.d0)
  nij(:) = - nij(:)
  tij(1) =  nij(2)*vv(n,3) - nij(3)*vv(n,2)
  tij(2) =  nij(3)*vv(n,1) - nij(1)*vv(n,3)
  tij(3) =  nij(1)*vv(n,2) - nij(2)*vv(n,1)
  tij(:) =  - tij(:) / dsqrt(tij(1)**2.d0+tij(2)**2.d0+tij(3)**2.d0)

  ee(1) =  (yp(mm(3))-yp(mm(1)))*(zp(mm(2))-zp(mm(1))) &
          -(zp(mm(3))-zp(mm(1)))*(yp(mm(2))-yp(mm(1)))
  ee(2) =  (zp(mm(3))-zp(mm(1)))*(xp(mm(2))-xp(mm(1))) &
          -(xp(mm(3))-xp(mm(1)))*(zp(mm(2))-zp(mm(1)))
  ee(3) =  (xp(mm(3))-xp(mm(1)))*(yp(mm(2))-yp(mm(1))) &
          -(yp(mm(3))-yp(mm(1)))*(xp(mm(2))-xp(mm(1)))
  area  = 0.5d0 * dsqrt(ee(1)**2.d0+ee(2)**2.d0+ee(3)**2.d0)
    
  ee(1) = (yp(3)-yp(4))*(zp(2)-zp(4)) - (zp(3)-zp(4))*(yp(2)-yp(4))
  ee(2) = (zp(3)-zp(4))*(xp(2)-xp(4)) - (xp(3)-xp(4))*(zp(2)-zp(4))
  ee(3) = (xp(3)-xp(4))*(yp(2)-yp(4)) - (yp(3)-yp(4))*(xp(2)-xp(4))
  area  = area + 0.5d0*dsqrt(ee(1)**2.d0+ee(2)**2.d0+ee(3)**2.d0)
  
  l = 2
  if(l==1) then
    ff(1) = ff(1) + (c * area - (wij- uij*area)*dot_product(nij,gv)*dtan(phi))
    ff(2) = ff(2) + wij*dot_product(tij,gv)
  end if
  if(l==2) then
    rb(1) = xx(i) ; rb(2) = yy(j) ; rb(3) = gl(i,j) - zvp(k)
    rb(:) = rb(:) - ctr(:)
    rg(1) = xx(i) ; rg(2) = yy(j) ; rg(3) = gl(i,j) - 0.5 * zvp(k)
    rg(:) = rg(:) - ctr(:)
    
    ee(1) = tij(2)*rb(3) - tij(3)*rb(2)
    ee(2) = tij(3)*rb(1) - tij(1)*rb(3)
    ee(3) = tij(1)*rb(2) - tij(2)*rb(1)
    a = dot_product(ee,vv(n,:))
    ff(1) = ff(1) + a * (c * area + (wij*(-1.d0*dot_product(nij,gv)) - uij*area)*dtan(phi))
    !ff(1) = ff(1) + a * (c * area + (wij*(-1.d0*dot_product(nij,gv)) &
    !        - uij*area*(dot_product(nij,gv)**2.d0))*dtan(phi))

    ee(1) = rb(2)*nij(3) - rb(3)*nij(2)
    ee(2) = rb(3)*nij(1) - rb(1)*nij(3)
    ee(3) = rb(1)*nij(2) - rb(2)*nij(1)
    a = dot_product(ee,vv(n,:))
    ee(1) = rg(2)*gv(3) - rg(3)*gv(2)
    ee(2) = rg(3)*gv(1) - rg(1)*gv(3)
    ee(3) = rg(1)*gv(2) - rg(2)*gv(1)
    b =  dot_product(ee,vv(n,:))
    ff(2) = ff(2) - wij*(a*(-1.d0)*dot_product(nij,gv) + b)
  end if

end do

if(k > 1) then
  fsp(1:n0) = ff(1)/ff(2)
  !write(*,*) ff(1), ff(2)
  if(fsp(1)<0.d0 .and. ff(1)<0.d0) fsp(1:n0) = 0.d0
  if(fsp(1)<0.d0 .and. ff(2)<0.d0) fsp(1:n0) = 10.d0
  if(fsp(1) > 1.d0) zvp(:) = 0.d0
end if

do k = 1,n0
  i = nep(k,1) ; j = nep(k,2)
  fs(i,j) = min(fs(i,j),fsp(k))
  zvol(i,j) = max(zvol(i,j),zvp(k))
end do
!!-------------------------------------------------------------------


  end do
end do
!$OMP end do
!$OMP end parallel



fname1 = './output/sfactor_'
write(fname,'(a,i4.4,a4)') trim(adjustl(fname1)), nf, '.txt'
open(11,file=fname,status='replace')
!fname1 = './output/zvol_'
!write(fname,'(a,i4.4,a4)') trim(adjustl(fname1)), nf, '.txt'
!open(12,file=fname,status='replace')
do i = 1,nx-1
  write(11,*) (real(fs(i,j)), j=1,ny-1)
  !write(12,*) (real(zvol(i,j)), j=1,ny-1)
end do
close(11)
!close(12)
end subroutine
! - - - - - - - - - - - - -


