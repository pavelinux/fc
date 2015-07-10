module variables 
 implicit none   
 integer :: nat,nn,npasos
 real*8  :: rho,lx,ly,dx,dy,upot,ukin,utot,temp
 real*8  :: sgm=1.0,eps=1.0,rcut,delta_t
 real*8,allocatable,dimension(:)::rx,ry
 real*8,allocatable,dimension(:)::vx,vy
 real*8,allocatable,dimension(:)::fx,fy
end module variables

!==================================
program md
use variables
implicit none
call datos  
call salva  
call fza    
call mdloop 
end program md

!====================================

subroutine datos
use variables
implicit none
integer:: i,j,k
real*8 :: rnd,mom_x,mom_y

open(1,file='run.txt',status='old',action='read')
 read(1,*)npasos
 read(1,*)nat
 read(1,*)rho
 read(1,*)delta_t
 read(1,*)rcut
 read(1,*)temp
close(1)

allocate(rx(nat),ry(nat))     
allocate(vx(nat),vy(nat))       
allocate(fx(nat),fy(nat))       
lx  = sqrt(dble(nat)/rho)       
ly  = lx                         
nn = sqrt(dble(nat))+1           
dx = lx/dble(nn)
dy = ly/dble(nn)


k=0
do i=1,nn
 do j=1,nn
   if(k<nat)then
     k=k+1
     rx(k) = dble(i)*dx
     ry(k) = dble(j)*dy
   end if
 end do
end do

mom_x =0.0
mom_y =0.0
!======================
do i=1,nat
 call random_number(rnd)
 vx(i) = 2.0*rnd - 1.0
 call random_number(rnd)
 vy(i) = 2.0*rnd - 1.0
 mom_x=mom_x+vx(i)
 mom_y=mom_y+vy(i)
end do
!========================
!stop
mom_x=mom_x/dble(nat)
mom_y=mom_y/dble(nat)
do i=1,nat
 vx(i) = vx(i) - mom_x 
 vy(i) = vy(i) - mom_y
end do
end subroutine datos

!=================================

subroutine salva
use variables
implicit none

integer :: i

OPEN(1,file='posi.xyz', status='unknown',action='write')

 write(1,*)nat
 write(1,*)'comentario'
 do i=1,nat
  write(1,*)'C',rx(i),ry(i),0.0
end do
!close(1)
end subroutine salva

!=====================================

subroutine fza
use variables
implicit none

integer ::i,j
real*8 ::rij,dulj,ulj,ucut

do i=1,nat
 fx(i) = 0.0 
 fy(i) = 0.0
end do
ucut = 4.0*eps*((sgm/rcut)**12 - (sgm/rcut)**6)!lennard-Jones

upot = 0.0 
do i=1,nat-1
 do j=i+1,nat
  dx = rx(i) - rx(j)
  dy = ry(i) - ry(j)

!=======condiciones==========================

 if(dx>0.50*lx)then
  dx = dx-lx
 elseif(dx<-0.50*lx)then
  dx = dx+lx
 end if
 if(dy>0.50*ly)then
  dy = dy-ly
 elseif(dy<-0.50*ly)then
  dy = dy+ly
 end if
  
  rij = sqrt(dx*dx+dy*dy)
  if(rij<=rcut)then
    ulj  = 4.0*eps*(sgm/rij)**6*((sgm/rij)**6-1.0)-ucut
    upot = upot + ulj
    dulj = 48.0*eps*(sgm/rij)**6*((sgm/rij)**6-0.50)/(rij*rij)
    fx(i) = fx(i) + dulj*dx
    fy(i) = fy(i) + dulj*dy
    fx(j) = fx(j) - dulj*dx
    fy(j) = fy(j) - dulj*dy
  end if
 end do
end do

end subroutine fza

!=====================================

subroutine mdloop
use variables

implicit none
integer,dimension(:),allocatable::frec
real*8, dimension(:),allocatable::g
integer ::i,paso,cont
real*8::fac,tins
integer::bin,no_bins,j
real*8::box_m,pi,da,delta,rij,lsup,linf

call salva_vel(0)
pi      = acos(-1.0)
delta   = 0.02
box_m   = min(lx*0.50,ly*0.50)
no_bins = int(box_m/delta)
allocate(frec(0:no_bins-1),g(0:no_bins-1))
frec=0
cont=0

do paso=1,npasos
 do i=1,nat
  vx(i)  =  vx(i) + 0.5*delta_t*fx(i)
  vy(i)  =  vy(i) + 0.5*delta_t*fy(i)
  rx(i)  =  rx(i) +     delta_t*vx(i)
  ry(i)  =  ry(i) +     delta_t*vy(i)
  if(rx(i)<0.0) rx(i) = rx(i) + lx
  if(ry(i)<0.0) ry(i) = ry(i) + ly
  if(rx(i)>lx ) rx(i) = rx(i) - lx
  if(ry(i)>ly ) ry(i) = ry(i) - ly
 end do
 call fza
 ukin=0.0
 do i=1,nat
  vx(i)  =  vx(i) + 0.5*delta_t*fx(i)
  vy(i)  =  vy(i) + 0.5*delta_t*fy(i)
  ukin   =  ukin  + vx(i)**2 + vy(i)**2
 end do
!---------------------------------------
 if(mod(paso,10)==0)then!mod(A.B) DICE A/B SI EL RESIDUO ES CERO THEN
  cont = cont+1
  do i=1,nat-1
   do j=i+1,nat
    dx = rx(i)-rx(j)
    dy = ry(i)-ry(j)
    if(dx>0.50*lx)then
     dx = dx-lx
    elseif(dx<-0.50*lx)then
     dx = dx+lx
    end if
    if(dy>0.50*ly)then
     dy = dy-ly
    elseif(dy<-0.50*ly)then
     dy = dy+ly
    end if
    rij = sqrt(dx*dx+dy*dy)
    bin = int(rij/delta)
    if(bin<no_bins)then
     frec(bin)=frec(bin)+2
    endif
   enddo
  enddo
 endif
!--------------------------------------------
 if(mod(paso,npasos/2)==0) call salva_vel(1)
 if(mod(paso,100)==0) call salva
 ukin = ukin*0.50  
 utot = upot + ukin
 tins = ukin/dble(nat)!temp instantanea
 fac  = sqrt(temp/tins)
 do i=1,nat
  vx(i)=vx(i)*fac
  vy(i)=vy(i)*fac
 end do
! write(*,'(i7,4f12.6)')paso,upot/dble(nat),ukin/dble(nat),utot/dble(nat),tins
 write(2,'(i7,4f12.6)')paso,upot/dble(nat),ukin/dble(nat),utot/dble(nat),tins
!=====================================
enddo

open (41,file="rdf.dat",action='write')
 do i=0,no_bins-1
  linf = dble(i)*delta
  lsup = linf+delta
  da   = pi*(lsup**2-linf**2)!PI*(r+dr)^2-pi*r^2
  g(i) = (lx*ly*frec(i))/(da*dble(nat**2))
  write(41,*)i*delta,g(i)/dble(cont)
 enddo
close(41)

 call salva_vel(2)
end subroutine mdloop

!=================================================
subroutine salva_vel(flag)
use variables
implicit none
integer::flag,i

 if(flag==0)then
  open(11,file='vel_ini.dat')
  else if (flag==1)then
  open(11,file='vel_med.dat')
  else
  open(11,file='vel_fin.dat')
 end if
do i=1,nat
 write(11,*)vx(i),vy(i)
enddo
close(11)
end subroutine salva_vel
