      PROGRAM LEE_ARCH
      IMPLICIT NONE
      INTEGER :: I, J
      REAL*8 :: A(4), B(4), C(4,2)

!      OPEN(UNIT=1, FILE='datos.dat', STATUS='OLD', ACTION='READ')
!      READ(1,*)
!      DO I = 1, 4
!        READ(1,*) A(I), B(I)
!        WRITE(*,*) A(I) + B(I)
!      END DO
!      CLOSE(1)

      OPEN(UNIT=1, FILE='datos.dat', STATUS='OLD', ACTION='READ')
      READ(1,*)
      DO I = 1, 5
        READ(1,*) (C(I,J), J=1,2)
        WRITE(*,*) (C(I,J), J=1,2)
      END DO
      
      CLOSE(1)

      END PROGRAM LEE_ARCH
