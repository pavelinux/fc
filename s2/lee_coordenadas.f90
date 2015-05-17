      PROGRAM LEE_ARCH
      IMPLICIT NONE
      INTEGER :: I, J
      REAL*8 :: A(1), B(1), C(1,2)
      OPEN(UNIT=1, FILE='datos.dat', STATUS='OLD', ACTION='READ')
      READ(1,*)
      DO I = 1, 5
        READ(1,*) (C(I,J), J=1,2)
        WRITE(*,*) (C(I,J), J=1,2)
      END DO
      
      CLOSE(1)

      END PROGRAM LEE_ARCH
