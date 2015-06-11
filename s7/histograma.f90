PROGRAM HISTOGRAMA
IMPLICIT NONE
REAL*8 :: DELTA,LI,LS,VX
INTEGER, DIMENSION(:), ALLOCATABLE::FREQ
INTEGER :: I, BIN, NO_BINS, NAT=1000
LS = 2.0
LI = 0.0
NO_BINS = 100
DELTA = (LS - LI) / NO_BINS
ALLOCATE(FREQ(NAT))
DO I = 1, NAT
    OPEN(1,FILE='datos_histograma.txt',STATUS='OLD',ACTION='READ')
    READ(1, *) VX
    VX = VX + 1.0
    BIN = VX / DELTA
    FREQ(BIN) = FREQ(BIN) + 1
END DO
CLOSE(1)

DO I = 0, NO_BINS
    OPEN(2,FILE='salida_histograma.dat',STATUS='UNKNOWN',ACTION='WRITE')
    WRITE(2,*) I * DELTA - 1, FREQ(I)
END DO
CLOSE(2)

END PROGRAM HISTOGRAMA
