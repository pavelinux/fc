MODULE VARIABLES
IMPLICIT NONE
INTEGER :: N,NN
REAL*8 RHO,LX,LY,DX,DY 
REAL*8, DIMENSION(:), ALLOCATABLE::RX,RY
END MODULE VARIABLES

PROGRAM MD
USE VARIABLES
IMPLICIT NONE
CALL COORDENADAS
CALL CELDA
END PROGRAM MD

SUBROUTINE COORDENADAS
USE VARIABLES
IMPLICIT NONE
INTEGER :: I,J,K
WRITE(*,*) 'DAME EL NUMERO DE ATOMOS'
READ(*,*) N
ALLOCATE(RX(N),RY(N))
WRITE(*,*) 'DAME LA DENSIDAD'
READ(*,*) RHO
LX = SQRT(DBLE(N)/RHO)
LY = LX
NN = SQRT(DBLE(N)) + 1
DX = LX/DBLE(NN)
DY = LY/DBLE(NN)
K=0
DO I=0, N
 DO J=0,N
 K=K+1
 IF(K<=N) THEN
     RX(K) = DBLE(J) * DX
     RY(K) = DBLE(J) * DY
 END IF
 END DO
END DO
END SUBROUTINE COORDENADAS

SUBROUTINE CELDA
USE VARIABLES
IMPLICIT NONE
INTEGER :: I
OPEN(1,FILE='COORD.xyz', STATUS='UNKNOWN', ACTION='WRITE')
WRITE(1,*) N
WRITE(1,*)
DO I=1,N
    WRITE(1,*)'C ', RX(I), RY(I), 0.0
END DO
CLOSE(1)
END SUBROUTINE CELDA
