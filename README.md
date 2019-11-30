# knnring
## k-NN algorithm using MPI in synchronous and asynchronous data tranfer

* To make MPI executable file *main_mpi* using *synchronous* method navigate to this directory and type: **make synchronous**.
* To make MPI executable file *main_mpi* using *asynchronous* method navigate to this directory and type: **make asynchronous**.
This also generates sequential executable *main_sequential*.

* Run using command: **mpirun -np 4 main_mpi**
