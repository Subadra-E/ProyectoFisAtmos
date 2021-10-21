
#include "timeIntegrator.hpp"




TimeIntegrator::TimeIntegrator()
{

}

//! Advance the state vector in time one time step
void TimeIntegrator::integrate( Euler1D &ps, double dx, double dt )
{
  // Runge-Kutta 2nd order

  // create two vector width q.getDim() components
  mVector<gridFunc1D> k1( ps.stateVector().getDim() );
  mVector<gridFunc1D> k2( ps.stateVector().getDim() );

  // give memory to each vector component
  for( int j=0; j<ps.stateVector().getDim(); j++ ){
    k1[j].create( ps.stateVector()[0].Npoints() );
    k2[j].create( ps.stateVector()[0].Npoints() );
  }

  k1 = dt * ps.RHS( ps.stateVector(), dx );
  boundaryCondition( k1 );
  k2 = dt * ps.RHS( ps.stateVector() + 0.5*k1, dx );
  boundaryCondition( k2 );

  ps.stateVector() = ps.stateVector() + k2;
  boundaryCondition( ps.stateVector() );
}
