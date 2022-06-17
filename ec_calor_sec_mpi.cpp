#include <mpi.h>
#include <math.h>
#include <ctime>
#include <cstdlib>  
#include <iostream>

double inicial ( double x, double tiempo ) {
  double limite;
  limite = 95.0;
  return limite;
}

double frontera ( double x, double tiempo ){
  double limite;
  if ( x < 0.5 ) limite = 100.0 + 10.0 * sin ( tiempo );
  else limite = 75.0;
  return limite;
}

int main(int argc, char **argv) {
  
  int rank{}, size{};
  int i,j,j_min = 0,j_max = 400,tag,n = 10;
  int p; //new parameters
  double k = 0.002;
  double tiempo,dt,tmax = 10.0,tmin = 0.0,tnew;
  double u[n],unew[n],x[n],dx;
  double x_max = 1.0,x_min = 0.0;
  dt = ( tmax - tmin ) / ( double ) ( j_max - j_min );
  dx = ( x_max - x_min ) / ( double ) ( n - 1 );
	x[0]=0;
  MPI_Init(&argc, &argv);                   // Initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank of the processor
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Total number of processors
  double btime{},etime{};
  p = size;
  if(rank == 0){
    std::cout<<"@∫dx Ecuacion de calor."<<std::endl;
    btime = MPI_Wtime();
  }
  
  for ( i = 1; i <= n + 1; i++ ){
  x[i] = ( ( double ) (rank * n + i - 1 ) * x_max +
  ( double ) ( p * n - rank * n - i) * x_min ) 
  / ( double ) ( p * n - 1 );
  }

  tiempo = tmin;
  u[0] = 0.0;
  for ( i = 1; i <= n; i++ ) u[i] = inicial ( x[i], tiempo );
  u[n+1] = 0.0;

  for ( j = 1; j <= j_max; j++ ){
    tnew += dt;
    if(rank > 0){
      tag = 201;
      MPI_Send( &u[1] , 1 , MPI_DOUBLE , rank-1 , tag , MPI_COMM_WORLD);  
    }
    if(rank < size-1){
      tag = 201;
      MPI_Recv( &u[n+1] , 1 , MPI_DOUBLE , rank+1 , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
    }
    if( rank < size -1){
      tag = 402;
      MPI_Send( &u[n] , 1 , MPI_DOUBLE , rank+1 , tag , MPI_COMM_WORLD);
    }
    if(rank > 0){
      tag = 402;
      MPI_Recv( &u[0] , 1 , MPI_DOUBLE , rank-1 , tag , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
    }

    for ( i = 1; i <= n; i++ ){
      unew[i] = u[i] 
      + ( dt * k / dx / dx ) * ( u[i-1] - 2.0 * u[i] + u[i+1] ); 
    }

    // putting the boundary values according to the boundary conditions.
    if(rank == 0) unew[1] = frontera( x[1], tnew );
	  if(rank == size-1) unew[n] = frontera( x[n], tnew );
    tiempo = tnew;

    for ( i = 1; i <= n; i++ ){
      	u[i] = unew[i];
	      if(j==j_max)printf(" %f %f\n",x[i],u[i]);
    }
  }

  MPI_Finalize();
}