program suma_mpi
implicit none
include "mpif.h"

integer :: n = 100, i, suma

integer :: nprocs, master = 0, esclavo, error,inicio, fin

call mpi_init(error)
call mpi_comm_size(mpi_comm_world,nprocs, error)
call mpi_comm_rank(mpi_comm_world,esclavo, error)
inicio = (esclavo * n) / nprocs + 1
fin = (esclavo + 1) * n / nprocs
write(*,*) inicio, fin
do i = inicio, fin
	suma = suma + 1
end do
!if(esclavo == master) write(*,*) suma
call mpi_finalize(error)
end program suma_mpi
