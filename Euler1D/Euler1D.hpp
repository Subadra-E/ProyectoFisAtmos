#ifndef Euler1D_H
#define Euler1D_H

#include "parameterReader.hpp"
#include "mvector.hpp"

class Euler1D
{
private:
  mVector<gridFunc1D> q;
  mVector<gridFunc1D> Flux;
  mVector<gridFunc1D> r_m, r_0, r_p, alpha, lambda, fluxRHS;
  gridFunc1D moment;
  gridFunc1D energ;
  int convFactor;

public:
  Euler1D( const parameterReader &p );
  mVector<gridFunc1D> &stateVector();
  const gridFunc1D &density() const;
  const gridFunc1D &vel() const;
  const gridFunc1D &pres() const;
  const gridFunc1D &momentum();
  const gridFunc1D &energy();
  void initialData( gridFunc1D &x, const parameterReader &p );
  mVector<gridFunc1D> RHS( const mVector<gridFunc1D> &iq, double dx );
  void output( gridFunc1D &x, double t );
};

inline double energy_PW( double d, double v, double p );
inline double momentum_PW( double d, double v );
void convergenceOutput( gridFunc1D &x, gridFunc1D &dens, gridFunc1D &xc, gridFunc1D &densc, int convFactor );
double minmod( double a, double b );
double superbee( double theta );
//void calc_dq( int i, int k, gridFunc1D &den, gridFunc1D &mom, gridFunc1D &ene, double *alpha, double *r_m, double *r_0, double *r_p, double *lambda );
void boundaryCondition( mVector<gridFunc1D> & q );
void calc_dq( int i, const mVector<gridFunc1D> &q, 
	      mVector<gridFunc1D> &alpha, 
	      mVector<gridFunc1D> &r_m, 
	      mVector<gridFunc1D> &r_0, 
	      mVector<gridFunc1D> &r_p,
	      mVector<gridFunc1D> &lambda );
void calcFlux3( const mVector<gridFunc1D> &iq, 
		mVector<gridFunc1D> &alpha,
		mVector<gridFunc1D> &r_m,
		mVector<gridFunc1D> &r_0,
		mVector<gridFunc1D> &r_p,
		mVector<gridFunc1D> &lambda,
		mVector<gridFunc1D> &Flux );

#endif
