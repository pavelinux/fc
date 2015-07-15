MODULE VARIABLES
IMPLICIT NONE
INTEGER :: N,NN,NPASOS
REAL*8 :: RHO,LX,LY,DX,DY,RIJ,M_X,M_Y,DELTA_T,R_CUT,TEMP
REAL*8 :: UPOT,UKIN,UTOT
REAL*8 :: ULJ,SGM = 1.0, EPS=1.0
REAL*8, DIMENSION(:), ALLOCATABLE::RX,RY
REAL*8, DIMENSION(:), ALLOCATABLE::VX,VY
REAL*8, DIMENSION(:), ALLOCATABLE::FX,FY
REAL*8  RND
END MODULE VARIABLES
!==================================
PROGRAM MD_DP
USE VARIABLES
IMPLICIT NONE
CALL COORDENADAS
CALL CELDA
CALL FUERZAS
CALL MDLOOP
END PROGRAM MD_DP
!==================================

SUBROUTINE COORDENADAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I,J,K

OPEN(1,FILE='input-p.txt',STATUS='OLD',ACTION='READ')
READ(1,*)
READ(1,*)
READ(1,*)NPASOS
READ(1,*)N
READ(1,*)RHO
READ(1,*)DELTA_T
READ(1,*)R_CUT
READ(1,*)TEMP
READ(1,*)LY
CLOSE(1)

ALLOCATE(RX(N),RY(N))
ALLOCATE(VX(N),VY(N))
ALLOCATE(FX(N),FY(N))

LX = DBLE(N) / (LY * RHO)
NN = SQRT(DBLE(N)) + 1 
DX = LX/DBLE(NN)
DY = LY/DBLE(NN)

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
!==================================
SUBROUTINE CELDA
USE VARIABLES
IMPLICIT NONE
INTEGER :: I
OPEN(2,FILE='coordenadas.xyz', STATUS='UNKNOWN', ACTION='WRITE')
WRITE(2,*) N, LX, LY
WRITE(2,*)
DO I=1,N
    WRITE(2,*)'C ', RX(I), RY(I), 0.0
END DO
!CLOSE(2) ! <- Descomentar para generar frames
END SUBROUTINE CELDA
!==================================
SUBROUTINE FUERZAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, J
REAL*8 :: DULJ, UCUT

DO I=1, N
    FX(I) = 0.0
    FY(I) = 0.0
END DO

UCUT = 4.0 * EPS * ((SGM/R_CUT)**12 - (SGM/R_CUT)**6) 
UPOT = 0.0

DO I = 1, N - 1
 DO J = I + 1, N 
     DX = RX(I) - RX(J)
     DY = RY(I) - RY(J)
    ! MINIMA IMAGEN
    IF(DX > 0.50 * LX) THEN
        DX = DX - LX
    ELSEIF(DX < - 0.50 * LX) THEN
        DX = DX + LX
    ENDIF 

    IF(DY > 0.50 * LY) THEN
        DY = DY - LY
    ELSEIF(DY < - 0.50 * LY) THEN
        DY = DY + LY
    ENDIF 

     RIJ = SQRT(DX * DX + DY * DY)
      IF (RIJ<= R_CUT) THEN
           ULJ = 4.0 * EPS * (SGM/RIJ)**6 * ((SGM/RIJ)**6 - 1.0) - UCUT
           UPOT = UPOT + ULJ
           DULJ = 48.0 * EPS * (SGM/RIJ)**6 * ((SGM/RIJ)**6 - 0.50)/(RIJ * RIJ)
           FX(I) = FX(I) + DULJ * DX
           FY(I) = FY(I) + DULJ * DY
           FX(J) = FX(J) - DULJ * DX
           FY(J) = FY(J) - DULJ * DY
        END IF
 END DO 
END DO
END SUBROUTINE FUERZAS
!==================================
SUBROUTINE MDLOOP
USE VARIABLES
IMPLICIT NONE
INTEGER :: I, PASO,CONT
INTEGER, DIMENSION(:), ALLOCATABLE::FREC
INTEGER,DIMENSION(:), ALLOCATABLE::NPARTX
INTEGER,DIMENSION(:), ALLOCATABLE::NPARTY
REAL*8, DIMENSION(:), ALLOCATABLE::G
REAL*8 :: TINS,FAC
INTEGER:: BIN,NO_BINS,J,CONT2=0
INTEGER:: MBINSX, MBINSY, BINX, BINY
REAL*8 :: BOX_M,PI,DA,DELTA,LSUP,LINF

PI = ACOS(-1.0)
DELTA = 0.02
BOX_M = MIN(LX * 0.50, LY * 0.50)
NO_BINS = INT(BOX_M / DELTA)
MBINSX = INT(LX / DELTA)
MBINSY = INT(LY / DELTA)
ALLOCATE(FREC(0:NO_BINS - 1),G(0:NO_BINS - 1))
ALLOCATE(NPARTX(0:MBINSX),NPARTY(0:MBINSY))
FREC = 0
CONT = 0

OPEN(3,FILE='energias.dat', STATUS='UNKNOWN', ACTION='WRITE') 
DO PASO = 1, NPASOS
    DO I= 1, N
        VX(I) = VX(I) + 0.50 * DELTA_T * FX(I)
        VY(I) = VY(I) + 0.50 * DELTA_T * FY(I)
        RX(I) = RX(I) + DELTA_T * VX(I)
        RY(I) = RY(I) + DELTA_T * VY(I)
        !BOUNDARY CONDITIONS
        IF(RX(I) < 0.0) RX(I) = RX(I) + LX
        IF(RY(I) < 0.0) RY(I) = RY(I) + LY
        IF(RX(I) > lx) RX(I) = RX(I) - LX
        IF(RY(I) > ly) RY(I) = RY(I) - LY
    END DO

    CALL FUERZAS
    UKIN=0.0

    DO I =1, N
        VX(I) = VX(I) + 0.50 * DELTA_T * FX(I)
        VY(I) = VY(I) + 0.50 * DELTA_T * FY(I)
        UKIN = UKIN + VX(I) ** 2 + VY(I)**2 
    END DO
    ! histograma perfil de rho_x, rho_y
    IF(MOD(PASO, 100) == 0) THEN
        CONT2 = CONT2 + 1
        DO I = 1, N
            BINX=RX(I) / DELTA
            BINY=RY(I) / DELTA
            IF(BINX <= MBINSX) NPARTX(BINX) = NPARTX(BINX) + 1
            IF(BINY <= MBINSY) NPARTY(BINY) = NPARTY(BINY) + 1
        END DO
    END IF
    ! histograma g(r)
    IF(MOD(PASO,10) == 0) THEN
        CONT = CONT + 1
        DO I = 1, N -1
            DO J = I + 1, N
                DX = RX(I)-RX(J)
                DY = RY(I)-RY(J)
                
                IF(DX > 0.50 * LX)THEN
                    DX = DX - LX
                ELSEIF(DX < -0.50 * LX) THEN
                    DX = DX + LX
                END IF
                IF(DY > 0.5 * LY)THEN
                    DY = DY - LY
                ELSEIF(DY < -0.50 * LY) THEN
                    DY = DY + LY
                END IF

                RIJ = SQRT(DX * DX + DY * DY)
                BIN = INT(RIJ / DELTA)

                IF(BIN < NO_BINS) THEN
                    FREC(BIN) = FREC(BIN) + 2
                END IF
            END DO
        END DO
    END IF    

    IF(MOD(PASO,10) == 0) CALL CELDA
    UKIN = UKIN * 0.50
    UTOT = UPOT + UKIN
    !Termostato
    TINS = UKIN / DBLE(N)
    FAC = SQRT(TEMP/TINS)

    DO I=1, N
    VX(I) = VX(I) * FAC
    VY(I) = VY(I) * FAC
    END DO

    WRITE(3,'(I7,X,4F12.6)')PASO,UPOT/DBLE(N),UKIN/DBLE(N),UTOT/(DBLE(N)),TINS
END DO
!================================== DO de la DINAMICA
! Escribe funcion de distribucion radial
OPEN(4, FILE='g_(r).dat', STATUS='UNKNOWN',ACTION='WRITE')
DO I = 0, NO_BINS - 1
    LINF = DBLE(I) * DELTA
    LSUP = LINF + DELTA
    DA = PI * ((LSUP)**2 - (LINF)**2)
    G(I) = LX*LY*FREC(I) / (DA * DBLE(N**2))
    WRITE(4,*) I* DELTA, G(I) / DBLE(CONT)
END DO
! ESCRIBE PERFILE DE DENSIDADES
OPEN(5, FILE='rho_x.dat',STATUS='UNKNOWN',ACTION='WRITE')
DO I=0, MBINSX
    WRITE(5,*) I*DELTA, NPARTX(I) / (DELTA * LY * CONT2)
END DO
CLOSE(5)
OPEN(6, FILE='rho_y.dat',STATUS='UNKNOWN',ACTION='WRITE')
DO I=0, MBINSY
    WRITE(6,*) I * DELTA, NPARTY(I) / (DELTA * LX * CONT2) 
END DO
CLOSE(6)
! CIERRA ARCHIVO DE G(R)
CLOSE(4)
CALL SALVA_VEL(2)
! cierra archivo de energias.dat
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
