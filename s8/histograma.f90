PROGRAM HISTOGRAMA
IMPLICIT NONE
REAL*8 :: DELTA,LI,LS,VX_I
INTEGER, DIMENSION(:), ALLOCATABLE::VX
INTEGER, DIMENSION(:), ALLOCATABLE::FREQ
INTEGER :: I, BIN, NO_BINS=100, NAT=1000, SUMA_FREQ

ALLOCATE(FREQ(NAT))
ALLOCATE(VX(NAT))

DO I = 1, NAT
    OPEN(1,FILE='vx_al_final.dat',STATUS='OLD',ACTION='READ')
    READ(1, *) VX(i)
END DO
CLOSE(1)

LS = MAXVAL(VX)
LI = MINVAL(VX)
DELTA = (LS - LI) / NO_BINS

DO I = 1, NAT
    OPEN(1,FILE='vx_al_final.dat',STATUS='OLD',ACTION='READ')
    READ(1,*) VX_I
    VX_i = VX_i + 1.0
    BIN = VX_i / DELTA
    write(*,*) vx_i, bin
    !FREQ(BIN) = FREQ(BIN) + 1
END DO
CLOSE(1)

suma_freq = 0.0
DO I = 0, NO_BINS
    OPEN(2,FILE='salida_histograma_vx_final.dat',STATUS='UNKNOWN',ACTION='WRITE')
    suma_freq = freq(i) + suma_freq
    WRITE(2,*) I * DELTA - 1, FREQ(I)
END DO
CLOSE(2)

if(suma_freq .ne. NAT) then
    write(*,*) 'wrong data. check your code!', suma_freq
else
    write(*,*) suma_freq
endif
END PROGRAM HISTOGRAMA
