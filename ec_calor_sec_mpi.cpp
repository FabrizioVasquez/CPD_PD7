#include <mpi.h>
#include <math.h>
#include <ctime>
#include <cstdlib>  
#include <iostream>


int main(int argc, char **argv) {
  
  int rank{}, size{};
  MPI_Init(&argc, &argv);                   // Initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank of the processor
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Total number of processors
  
  if(rank == 0){
    std::cout<<"@∫dx."<<std::endl;
  }

  MPI_Finalize();
}