PROGRAM SUMA_MPI
USE MPI
IMPLICIT NONE

INTEGER :: N = 100, I, SUMA
INTEGER :: status(mpi_comm_world) 
INTEGER :: NPROCS, MASTER = 0, ESCLAVO, ERROR,INICIO, FIN, resultado=0

CALL MPI_INIT(ERROR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, ERROR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, ESCLAVO, ERROR)

INICIO = (ESCLAVO * N) / NPROCS + 1
FIN = (ESCLAVO + 1) * N / NPROCS

DO I = INICIO, FIN
    SUMA = SUMA + i
END DO

if(esclavo == master) then
call mpi_send(suma, 1, mpi_integer, 1, 100, mpi_comm_world,error)
else
call mpi_recv(resultado, 1, mpi_integer,master,100,mpi_comm_world,status,error)
write(*,*) suma + resultado
end if


CALL MPI_FINALIZE(ERROR)

END PROGRAM SUMA_MPI
