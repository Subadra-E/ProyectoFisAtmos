#ifndef gridFunc1D_H
#define gridFunc1D_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;

enum boundaryType_t { PERIODIC, NEWMAN, DIRICHLET };

class gridFunc1D
{
private:
  int n_points;
  double *data;
  double *datamid;
  int boundaryType;
  // boundary type
  //-1 = no value has been assigned
  // 0 = periodic
  // 1 = newman (derivative = 0)
  // 2 = dirichlet
  int isFlux; // 0 or 1, for variables holding fluxes values at
  // semi integer grid points.
 

public:
  gridFunc1D();
  gridFunc1D( int );
  gridFunc1D( const gridFunc1D & );
  ~gridFunc1D();
  void create( int );
  void erase();
  double &operator[]( float ) const;
  double &operator[]( float );
  int Npoints() const;
  void setBoundaryCondition( int);
  void setIsFlux( int);
  void outputGnuplotFake( ofstream &, gridFunc1D &, const double);
  void outputGnuplot( ofstream &, gridFunc1D &, const double) const;
  void outputByLine( ofstream &, const double t ) const;  
  void outputByColumn( ofstream &, gridFunc1D &, const double t ) const;
  gridFunc1D operator+( const gridFunc1D &B ) const;
  gridFunc1D operator-( const gridFunc1D &B ) const;
  friend gridFunc1D operator-( const double a, const gridFunc1D &B );
  gridFunc1D operator*( const gridFunc1D &B ) const;
  gridFunc1D operator*( const double &b ) const;
  friend gridFunc1D operator*( const double &a, const gridFunc1D &B );
  gridFunc1D operator/( const double &b ) const;
  gridFunc1D operator/( const gridFunc1D &B ) const;
  friend gridFunc1D operator/( const double a, const gridFunc1D &B );
  const gridFunc1D &operator=( const gridFunc1D &B );
  const gridFunc1D &operator=( const double b );
  const gridFunc1D &operator=( const int b );
  friend gridFunc1D sqrt( const gridFunc1D &A );
  void ishuge( double a );
};

ostream &operator<<( ostream &os, const gridFunc1D &A );



#endif
