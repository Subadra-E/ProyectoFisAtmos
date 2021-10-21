
#include <fstream>
#include <iostream>
#include <string>

#include "gridFunc1D.hpp"
#include "Euler1D.hpp"
#include "parameterReader.hpp"
#include "timeIntegrator.hpp"

class Simulation
{
protected:
  const double fudge; //!< A small quantity
  int n_points_x; //!< Number of points in the x direction
  double xMin; //!< Minimum value of x coordinate
  double xMax; //!< Maximum value of x coordinate
  double CFL;  //!< Courant-Friederich-Levy factor
  double finalTime; //!< Final time, where simulation stops
  double outEveryTime; //!< How often in time we do output
  double time; //!< The time coordinate (variable)
  
  double dx; //!< Separation between points in x
  double dt; //!< Separation between points in time
  int outEvery; //!< How many iterations we do output
  int ITMAX; //!< Maximum nunmber of iterations, when simulation stops
  gridFunc1D x; //!< x coordinate
  Euler1D *phys; //!< The physical system of equations
  TimeIntegrator tInteg; //!< The time integration algorithm

 
public:
  Simulation();
  Simulation( const parameterReader &p );
  void initialData( const parameterReader &p );
  void evolve();
};

Simulation::Simulation() : fudge(1.0e-10), n_points_x(100), xMin(0.0), xMax(1.0), CFL(0.1), finalTime(10.0), outEveryTime(0.1), time(0.0)
{
}

//Simulation::Simulation( int Npx, double xmin, double xmax, double cfl, double finaltime, double outEveTime, int cf ) : fudge(1.0e-10)
Simulation::Simulation( const parameterReader &p ) : fudge(1.0e-10)
{
  // variables initialized with information given by user
  n_points_x   = atoi( p.getParam( "n_points" ) );
  xMin         = atof( p.getParam( "x_min" ) );
  xMax         = atof( p.getParam( "x_max" ) );
  finalTime    = atof( p.getParam( "final_time" ) );
  outEveryTime = atof( p.getParam( "out_every_time" ) );
  CFL          = atof( p.getParam( "CFL" ) );
  time = 0.0;
  

  // variables initialized from variables above
  dx = (xMax - xMin) / n_points_x;
  dt = CFL * dx;
  outEvery = (int)(outEveryTime/dt+fudge);
  ITMAX = (int)(finalTime/dt+fudge);

  // initialize x coordinate
  x.create( n_points_x );
  for( int i=1; i<=n_points_x; i++ )
    x[i] = xMin + dx*(i-1);
  // x[i] = xMin + dx*i - 0.5*dx;

  // create the physical system
  phys = new Euler1D( p );

  // Time integrator is initializaed with a reference to the function
  // that computes the RHS of physical system phys.
  // tInteg = new TimeIntegrator();
}


void Simulation::initialData( const parameterReader &p )
{
  phys->initialData( x, p );

  cout <<"it= "<<0<<"\t"<<"time= "<<time<<endl;
  phys->output( x, time );
}


void Simulation::evolve()
{
  // integrating cycle
  for( int j=1; j<=ITMAX; j++ ){
    //phys->RHS( dt, dx );
    tInteg.integrate( *phys, dx, dt );
    time += dt;

    if ( j%outEvery==0 ) {
      cout <<"it= "<<j<<"\t"<<"time= "<<time<<endl;
      phys->output( x, time );
    }
  }
}





//=========================================================================

int main()
{
  // read parameter file
  parameterReader p( "Parameters" );

  // create simulation
  Simulation sim( p );

  // initial data
  sim.initialData( p );

  // evolution
  sim.evolve();

  return 0;
}

