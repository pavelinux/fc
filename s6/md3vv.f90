MODULE VARIABLES
IMPLICIT NONE
INTEGER :: N,NN,npasos=10000
REAL*8 :: RHO,LX,LY,LZ,DX,DY,DZ,RIJ,UPOT,ULJ,M_X,M_Y,delta_t=0.001, ukin
REAL*8 :: SGM = 1.0, EPS=1.0, RCUT=2.50
REAL*8, DIMENSION(:), ALLOCATABLE::RX,RY,RZ
REAL*8, DIMENSION(:), ALLOCATABLE::VX,VY,VZ
REAL*8, DIMENSION(:), ALLOCATABLE::FX,FY,FZ
REAL*8 RND
END MODULE VARIABLES

PROGRAM MD_F
USE VARIABLES
IMPLICIT NONE
CALL COORDENADAS
CALL CELDA
CALL FUERZAS
call mdloop
END PROGRAM MD_F

SUBROUTINE COORDENADAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: H,I,J,K,M
WRITE(*,*) 'DAME EL NUMERO DE ATOMOS'
READ(*,*) N
ALLOCATE(RX(N),RY(N),RZ(N))
ALLOCATE(VX(N),VY(N),VZ(N))
ALLOCATE(FX(N),FY(N),FZ(N))
!WRITE(*,*) 'DAME LA DENSIDAD'
!READ(*,*) RHO
RHO = 0.1
LX = SQRT(DBLE(N/RHO))
LY = LX
LZ = LX
NN = SQRT(DBLE(N)) + 1 
DX = LX/DBLE(NN)
DY = DX
DZ = DX
K=0
DO H=1,NN 
    DO I=1,NN - 1
         DO J=1,NN - 1
             IF(K<N) THEN
                 K=K+1
                 RX(K) = DBLE(I) * DX
                 RY(K) = DBLE(J) * DY
                 RZ(K) = DBLE(H) * DZ
             END IF
         END DO
    END DO
END DO
M_X = 0.0
M_Y = 0.0
DO I =1, N
    CALL RANDOM_NUMBER(RND)
    VX(I) = 2.0 * RND - 1.0
    VY(I) = 2.0 * RND - 1.0
    M_X = M_X + VX(I)
    M_Y = M_Y + VY(I)
END DO
M_X = M_X / DBLE(N)
M_Y = M_Y / DBLE(N)
DO I = 1, N
    VX(I) = VX(I)-M_X
    VY(I) = VY(I)-M_Y
    VZ(I) = VZ(I)-M_Y
END DO
END SUBROUTINE COORDENADAS

SUBROUTINE CELDA
USE VARIABLES
IMPLICIT NONE
INTEGER :: I
OPEN(1,FILE='COORDENADAS.xyz', STATUS='UNKNOWN', ACTION='WRITE')
WRITE(1,*) N
WRITE(1,*)
DO I=1,N
    !WRITE(1,*)'C ', RX(I), RY(I), RZ(I)
    WRITE(1,*)'C ', RX(I), RY(I), 0.0
END DO
CLOSE(1)
END SUBROUTINE CELDA

SUBROUTINE FUERZAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, J
REAL*8 :: DULJ
DO I=1, N
    FX(I) = 0.0
    FY(I) = 0.0
!    FZ(I) = 0.0
END DO
UPOT = 0.0
DO I = 1, NN - 1
 DO J = I+1, NN
     DX = RX(I) - RX(J)
     DY = RY(I) - RY(J)
     RIJ = SQRT(DX * DX + DY * DY)
      IF (RIJ<= RCUT) THEN
           ULJ = 4.0 * EPS * ((SGM/RIJ)**6 * ( (SGM/RIJ)**6 - 1.0)) - UCUT
!           ULJ = 4.0 * EPS * ((SGM/RIJ)**6 * ( (SGM/RIJ)**6 - 1.0)) 
           UPOT = UPOT + ULJ
           DULJ = 48.0 * EPS * (SGM/RIJ)**6 * ((SGM/RIJ)**6 - 0.5)/(RIJ * RIJ)
           FX(I) = FX(I) + DULJ * DX
           FY(I) = FY(I) + DULJ * DY
           FX(J) = FX(J) - DULJ * DX
           FY(J) = Fy(J) - DULJ * DY
        END IF
 END DO 
END DO
!!WRITE(*,*) sum(fx),sum(y),ucut
WRITE(*,*) sum(fx),sum(fy)
stop
END SUBROUTINE FUERZAS

subroutine mdloop
use variables
implicit none
integer :: i, paso
do paso = 1, npasos
do i =1, N
vx(i) = vx(i) + 0.50 * delta_t * fx(i)
vy(i) = vy(i) + 0.50 * delta_t * fy(i)
!vz(i) = vz(i) + 0.50 * delta_t * fz(i)
rx(i) = rx(i) + delta_t *vx(i)
ry(i) = ry(i) + delta_t *vy(i)
!if(rx(i) < 0.0) rx(i) = rx(i) + lx
!if(ry(i) < 0.0) ry(i) = ry(i) + ly
!if(rx(i) > 0.0) rx(i) = rx(i) - lx
!if(ry(i) > 0.0) ry(i) = ry(i) - ly
!rz(i) = rz(i) + delta_t *vz(i)
end do
call fuerzas
ukin=0.0
do i =1, N
    vx(i) = vx(i) + 0.50 * delta_t * fx(i)
    vy(i) = vy(i) + 0.50 * delta_t * fy(i)
    !vz(i) = vz(i) + 0.50 * delta_t * fz(i)
    ukin = ukin + vx(i) ** 2 + vy(i)**2! + vz(i)**2
end do
if(mod(paso,100)== 0) call celda
    ukin = ukin * 0.50
    write(*,*)paso,upot/dble(N),ukin/dble(N)
    write(2,*)paso,upot/dble(N),ukin/dble(N)
end do

end subroutine mdloop
