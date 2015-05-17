      PROGRAM FIND_MIN
      IMPLICIT NONE
      INTEGER :: I
      REAL numero
      numero = 1.0
      DO I = 1, 10000000
          numero = numero / 2
          WRITE(*,*)  (numero), I
          IF (numero == 0.00) THEN
          STOP
          END IF
      END DO



      END PROGRAM FIND_MIN 
