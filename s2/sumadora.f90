      PROGRAM SUMADORA
     
      IMPLICIT NONE
      INTEGER*8 :: I, SUMA

      SUMA = 0
      DO I = 1, 500
        SUMA = I + SUMA
      END DO
      WRITE(*,*) "EL VALOR DE LA SUMA ES: ", SUMA

      END PROGRAM SUMADORA
