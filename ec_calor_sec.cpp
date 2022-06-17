# include <iostream>
# include <cmath>
# include <ctime>

double frontera ( double x, double tiempo ){
  double limite;
  if ( x < 0.5 ) limite = 100.0 + 10.0 * sin ( tiempo );
  else limite = 75.0;
  return limite;
}

double inicial ( double x, double tiempo ) {
  double limite;
  limite = 95.0;
  return limite;
}

int main ( int argc, char *argv[] )
{
  // j_max -> intervalos de tiempo
  // n     -> numero de intervalos
  // 100.0 + 10.0 *sin(tiempo) -> condicion de frontera inicial
  // limite = 75.0 -> limite extremo
  // k -> constante de multiplicacion
  // [x_min,x_max] -> dominio total
  // [t_min,t_max] -> intervalo total de tiempo
  // inicial() -> límite = 95.0 -> temperatura inicial


  int i,j,j_min = 0,j_max = 400,tag,n = 10;
  double k = 0.002;
  double tiempo,dt,tmax = 10.0,tmin = 0.0,tnew;
  double u[n],unew[n],x[n],dx;
  double x_max = 1.0,x_min = 0.0;
  dt = ( tmax - tmin ) / ( double ) ( j_max - j_min );
  dx = ( x_max - x_min ) / ( double ) ( n - 1 );
	x[0]=0;
	for ( i = 1; i < n; i++ ){
		x[i]=x[i-1]+dx;
  	}

//  Inicializacion.
  tiempo = tmin;
  u[0] = 0.0;
  for ( i = 1; i <= n; i++ ) u[i] = inicial ( x[i], tiempo );
  u[n+1] = 0.0;

// Valores de temperatura al siguiente intervalo de tiempo
  for ( j = 1; j <= j_max; j++ ){
    tnew+=dt; 

// actualizar temperatura
    for ( i = 1; i <= n; i++ ){
      unew[i] = u[i] 
      + ( dt * k / dx / dx ) * ( u[i-1] - 2.0 * u[i] + u[i+1] ); 
    }


// aplicar condiciones de frontera
	    unew[1] = frontera( x[1], tnew );
	    unew[n] = frontera( x[n], tnew );


// Actualizar tiempo y temperatura
    tiempo = tnew;

    for ( i = 1; i <= n; i++ ){
      	u[i] = unew[i];
	if(j==j_max)printf("%f %f %f\n",tiempo,x[i],u[i]);
    }

  }

  return 0;
}


