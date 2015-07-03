MODULE VARIABLES
IMPLICIT NONE
INTEGER :: N,NN,NPASOS
REAL*8 :: RHO,LX,LY,DX,DY,RIJ,UPOT,UKIN,ukin_t,UTOT,ULJ,M_X,M_Y,DELTA_T,R_CUT,TEMP
REAL*8 :: SGM = 1.0, EPS=1.0
REAL*8, DIMENSION(:), ALLOCATABLE::RX,RY
REAL*8, DIMENSION(:), ALLOCATABLE::VX,VY
REAL*8, DIMENSION(:), ALLOCATABLE::FX,FY
integer, dimension(:), allocatable::histo_bin
real*8, dimension(:), allocatable::g
REAL*8  RND
REAL*8 :: T_I,T_F
END MODULE VARIABLES

MODULE UNIDADES_REDUCIDAS
IMPLICIT NONE
! parametros lj para el Ar
REAL*8 :: S_AR = 3.4E-10, EPS_AR=1.65E-21,M_AR=6.69E-26
REAL*8 :: T_R, L_R, RHO_R, UTOT_R
! K-Boltzman
REAL*8 :: kB = 1.3806E-23
END MODULE UNIDADES_REDUCIDAS

PROGRAM MD_F
USE VARIABLES
IMPLICIT NONE
CALL CPU_TIME(T_I)
CALL COORDENADAS
CALL CELDA
CALL FUERZAS
CALL MDLOOP
CALL CPU_TIME(T_F)
!CALL LOGFILE
END PROGRAM MD_F

SUBROUTINE COORDENADAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I,J,K,M

OPEN(1,FILE='input.txt',STATUS='OLD',ACTION='READ')
READ(1,*)
READ(1,*)
READ(1,*)NPASOS
READ(1,*)N
READ(1,*)RHO
READ(1,*)DELTA_T
READ(1,*)R_CUT
READ(1,*)TEMP
READ(1,*)
CLOSE(1)

ALLOCATE(RX(N),RY(N))
ALLOCATE(VX(N),VY(N))
ALLOCATE(FX(N),FY(N))

LX = SQRT(DBLE(N/RHO))
LY = LX
NN = SQRT(DBLE(N)) + 1 
DX = LX/DBLE(NN)
DY = lY/dble(NN)
K=0
DO I=1,NN 
     DO J=1,NN
         IF(K<N) THEN
             K=K+1
             RX(K) = DBLE(I) * DX
             RY(K) = DBLE(J) * DY
         END IF
     END DO
END DO
M_X = 0.0
M_Y = 0.0
DO I =1, N
    CALL RANDOM_NUMBER(RND)
    VX(I) = 2.0 * RND - 1.0
    CALL RANDOM_NUMBER(RND)
    VY(I) = 2.0 * RND - 1.0
    M_X = M_X + VX(I)
    M_Y = M_Y + VY(I)
END DO
M_X = M_X / DBLE(N)
M_Y = M_Y / DBLE(N)
DO I = 1, N
    VX(I) = VX(I) - M_X
    VY(I) = VY(I) - M_Y
END DO
END SUBROUTINE COORDENADAS

SUBROUTINE CELDA
USE VARIABLES
IMPLICIT NONE
INTEGER :: I
OPEN(2,FILE='coordenadas.xyz', STATUS='UNKNOWN', ACTION='WRITE')
WRITE(2,*) N
WRITE(2,*)
DO I=1,N
    WRITE(2,*)'C ', RX(I), RY(I), 0.0
END DO
CLOSE(2)
END SUBROUTINE CELDA

SUBROUTINE FUERZAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, J
REAL*8 :: DULJ, ucut

DO I=1, N
    FX(I) = 0.0
    FY(I) = 0.0
END DO
UCUT = 4.0 * EPS * ((SGM/R_CUT)**12 - (SGM/R_CUT)**6) 
UPOT = 0.0
DO I = 1, NN - 1
 DO J = I + 1, N 
     DX = RX(I) - RX(J)
     DY = RY(I) - RY(J)
    ! MINIMA IMAGEN
    IF(DX > 0.5 * LX) THEN
        DX = DX - LX
    ELSEIF(DX < - 0.5 * LX) THEN
        DX = DX + LX
    ENDIF 
    IF(DY > 0.5 * LY) THEN
        DY = DY - LY
    ELSEIF(DY < - 0.5 * LY) THEN
        DY = DY + LY
    ENDIF 
    !
     RIJ = SQRT(DX * DX + DY * DY)
      IF (RIJ<= R_CUT) THEN
           ULJ = 4.0 * EPS * ((SGM/RIJ)**6 * ( (SGM/RIJ)**6 - 1.0)) - UCUT
           UPOT = UPOT + ULJ
           DULJ = 48.0 * EPS * (SGM/RIJ)**6 * ((SGM/RIJ)**6 - 0.5)/(RIJ * RIJ)
           FX(I) = FX(I) + DULJ * DX
           FY(I) = FY(I) + DULJ * DY
           FX(J) = FX(J) - DULJ * DX
           FY(J) = FY(J) - DULJ * DY
        END IF
 END DO 
END DO
END SUBROUTINE FUERZAS

SUBROUTINE MDLOOP
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, PASO
REAL*8 :: TINS,FAC
OPEN(3,FILE='energias.dat', STATUS='UNKNOWN', ACTION='WRITE') 
!call gder(0)
DO PASO = 1, NPASOS
    DO I =1, N
        VX(I) = VX(I) + 0.50 * DELTA_T * FX(I)
        VY(I) = VY(I) + 0.50 * DELTA_T * FY(I)
        RX(I) = RX(I) + DELTA_T *VX(I)
        RY(I) = RY(I) + DELTA_T *VY(I)
        !BOUNDARY CONDITIONS
        IF(RX(I) < 0.0) RX(I) = RX(I) + LX
        IF(RY(I) < 0.0) RY(I) = RY(I) + LY
        IF(RX(I) > 0.0) RX(I) = RX(I) - LX
        IF(RY(I) > 0.0) RY(I) = RY(I) - LY
    END DO

    CALL FUERZAS
    UKIN=0.0

    DO I =1, N
        VX(I) = VX(I) + 0.50 * DELTA_T * FX(I)
        VY(I) = VY(I) + 0.50 * DELTA_T * FY(I)
        UKIN = UKIN + VX(I) ** 2 + VY(I)**2 
    END DO

    IF(MOD(PASO,100)== 0) CALL CELDA
    UKIN = UKIN * 0.50
    UTOT = UPOT + UKIN
    UKIN_T = 0.5 * UKIN
    TINS = UKIN_T / DBLE(N)
    FAC = SQRT(TEMP/TINS)

    DO I=1, N
    VX(I) = VX(I) * FAC
    VY(I) = VY(I) * FAC
    END DO

    WRITE(3,'(I7,X,4F12.6)')PASO,UPOT/DBLE(N),UKIN/DBLE(N),UTOT/(DBLE(N)),TINS
END DO
CALL SALVA_VEL(2)
!call gder(1)
CLOSE(3)
END SUBROUTINE MDLOOP

SUBROUTINE SALVA_VEL(FLAG)
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, FLAG
IF(FLAG == 0) THEN
    OPEN(7,FILE='velocidades_al_inicio.dat',STATUS='UNKNOWN',ACTION='WRITE')
    DO I = 1, N
    WRITE(7,*) VX(I), VY(I)
    END DO
    CLOSE(7)
ELSEIF(FLAG == 1) THEN
    OPEN(8, FILE='velocidades_a_mitad.dat',STATUS='UNKNOWN',ACTION='WRITE')
    DO I = 1, N
    WRITE(8,*) VX(I), VY(I)
    END DO
    CLOSE(8)
ELSE
    OPEN(9,FILE='velocidades_al_final.dat',STATUS='UNKNOWN', ACTION='WRITE')
    DO I = 1, N
    WRITE(9,*) VX(I), VY(I)
    END DO
    CLOSE(9)
END IF
END SUBROUTINE SALVA_VEL

subroutine gder
use variables
implicit none
integer :: k, flag
integer :: i,j,bin,num_bins
real*8 :: box_m, linf, lsup,da,A 
real*8 :: delta_bin = 0.02, box_min
REAL, PARAMETER :: Pi = 3.1415927
if(flag == 0) then
! calcula bin
! MINIMA IMAGEN
IF(DX > 0.5 * LX) THEN
    DX = DX - LX
ELSEIF(DX < - 0.5 * LX) THEN
    DX = DX + LX
ENDIF 
IF(DY > 0.5 * LY) THEN
    DY = DY - LY
ELSEIF(DY < - 0.5 * LY) THEN
    DY = DY + LY
ENDIF 
!
RIJ = SQRT(DX ** 2 + DY** 2)
do i = 1, N - 1
    do j = i + 1, N
! calcula no_bins
    bin = rij / delta_bin
    box_min = min(lx/2, ly/2)
    num_bins = box_min / delta_bin
    allocate(histo_bin(0:num_bins - 1))
    histo_bin = 0
    if(bin <= num_bins) then
        histo_bin(bin) = histo_bin(bin) + 2
    end if
    end do
end do
else
allocate(g(num_bins - 1))
open(11, file='rdf.dat', status='unknown',action='write')
do i = 0, num_bins - 1
    linf = i * delta_bin
    lsup = linf + delta_bin
    da = PI * ((lsup)**2 - (linf)**2)
    g(i) = lx*ly*histo_bin(i) / (N**2*da)
    write(11,*) i* delta_bin, g(i)
end do
close(11)
end if
end subroutine gder

