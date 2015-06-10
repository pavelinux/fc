PROGRAM HISTOGRAMA
IMPLICIT NONE
REAL*8 :: DELTA,LI,LS,VX
INTEGER, DIMENSION(:), ALLOCATABLE::FREQ
INTEGER :: I, BIN, NO_BINS, NAT=1000
LS = 2.0
LI = 0.0
NO_BINS = 100
DELTA = (LS - LI) / NO_BINS
DO I = 1, NAT
    OPEN(1,file='datos_histograma.txt',status='old',action='read')
    READ(1, *) VX
    VX = VX + 1.0
    BIN = VX / DELTA
    FREQ(BIN) = FREQ(BIN) + 1
    close(1)
END DO

DO I = 1, NO_BINS
    open(2,file='salida_histograma.dat',status='unknown',action='write')
    WRITE(2,*) I * DELTA - 1, FREQ(I)
END DO
close(2)

END PROGRAM HISTOGRAMA
