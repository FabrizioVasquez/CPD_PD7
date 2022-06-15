make:
	mpic++ ec_calor_sec_mpi.cpp
	mpirun -np 4 ./a.out
	rm ./a.out