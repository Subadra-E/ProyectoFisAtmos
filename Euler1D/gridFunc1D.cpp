
#include <iomanip>
#include "gridFunc1D.hpp"

//! Constructor with no arguments
//!
//! It initializes everything to zero
gridFunc1D::gridFunc1D()
{
  n_points = 0;
  data = NULL;
  datamid = NULL;
  boundaryType = -1;
  isFlux = 0;
}

//! Initializes and creates space to hold n elements
gridFunc1D::gridFunc1D( int n )
{
  n_points = 0;
  data = NULL;
  datamid = NULL;
  boundaryType = -1;
  isFlux = 0;

  create( n );
}

//! Copy constructor
gridFunc1D::gridFunc1D( const gridFunc1D &A )
{
  n_points = 0;
  data = NULL;
  datamid = NULL;
  boundaryType = -1;
  isFlux = 0;

  create( A.n_points );

  for( int i=0; i<n_points; i++ )
    data[i] = A.data[i];
}


gridFunc1D::~gridFunc1D()
{
  erase();
}

//! Creates and resize the objects
//!
//! This function assigns space and also can resize an already
//! existing object
//!
void gridFunc1D::create( int n )
{
  erase();

  if ( n > 0 ){
    n_points = n;
    data = new double[n];
    datamid = new double[n+1];

    for( int i=0; i<n; i++ )   data[i] = 0.0;
    for( int i=0; i<n+1; i++ ) datamid[i] = 0.0;
  }
  else {
    cout << "ERROR: the number of points must be positive."<<endl;
    exit(1);
  }
}

void gridFunc1D::erase()
{
  if ( data != NULL ) delete [] data;
  if ( datamid != NULL ) delete [] datamid;

  n_points = 0;
  data = NULL;
  datamid = NULL;
}

void gridFunc1D::setBoundaryCondition( int bc)
{
  if ( (bc == 0) || ( bc ==1 ) || ( bc == 2 ) ) 
    boundaryType = bc;
  else {
    cout << "ERROR: invalid boundary condition" << endl;
    exit(1);
  }
    
}



void gridFunc1D::setIsFlux( int a )
{
  if ((a==0) || (a==1))
    isFlux = a;
  else {
    cout << "ERROR: isFlux at of range" << endl;
    exit(1);
  }
}

//! Gnuplot style output (fake)
//!
//! This function is trick to make a 1D variable look like a 2D one.
//! Its purpose is the make a density plot with Gnuplot
//!

void gridFunc1D::outputGnuplotFake( ofstream &out, gridFunc1D &x, const double t ) 
{
  out << "#time = " << t << endl;

  for( int i=1; i<=n_points; i++ )
    out << 0.0 << "\t" << x[i] << "\t" << (*this)[i] << endl;

  out << endl;

  for( int i=1; i<=n_points; i++ )
    out << 1.0 << "\t" << x[i] << "\t" << (*this)[i] << endl;

  out << endl << endl;
}


//! Gnuplot style output
//!
//! First column is time, second column is spatial variable
//! and third column is the field itself
//!

void gridFunc1D::outputGnuplot( ofstream &out, gridFunc1D &x, const double t ) const
{
  for( int i=1; i<=n_points; i++ )
    out << t << "\t" << x[i] << "\t" << (*this)[i] << endl;

  out << endl << endl;
}




//! Output all the values of the variable in a single line
//!
//! Time is in the first column followed by all the pointwise
//! values of the variable. The position is not output. This is
//! useful the integrate the velocity to get the position of
//! hypothetical particles.
//!
void gridFunc1D::outputByLine( ofstream &out, const double t ) const
{
  out << t;

  for( int i=1; i<=n_points; i++ )
    out << "\t" << (*this)[i];

  out << endl;
}


//! ygraph output style
//!
//! ygraph output style consists in  starting a block with the current time
//! followed with the position and the value of the variable in one line
//!
void gridFunc1D::outputByColumn( ofstream &out, gridFunc1D &x, const double t ) const
{
  out << "#time = " << t << endl;

  for( int i=1; i<=n_points; i++ )
    out << x[i] << "\t" << std::setprecision(15) << (*this)[i] << endl;

  out << endl;
}




// we numerate the points from 1 the N
// This functions make the asumption that the grid function is
// periodic. This means that it accepts values of the index that
// are higher than N and lower than 1. For values grater than N
// it wraps around the values as a true periodic function. The
// same goes for values of the index lower than one.
double &gridFunc1D::operator[]( float j ) const
{
  int i;
  float checkInt = j - floorf(j);

  // check if j is integer
  if (checkInt == 0) {
    i=(int)j;
    if ( (1<=i) && (i<=n_points) )
      return data[i-1];
    // for i higher than N
    else if ( (n_points<i) && (i<2*n_points) && (boundaryType == 0) )
      return data[i-n_points-1];
    // for i lower than 1
    else if ( (-n_points<i) && (i<1) && (boundaryType == 0) )
      return data[n_points-i-1];
    else {
      cout << "ERROR: index out range"<<endl;
      cout << "Index requested = "<<i<<endl;
      exit(1);
    }
  }

  // check if j is not integer
  else if ( checkInt == 0.5 ){
    // check if we are at the boundaries
    if ((j==0.5) && (boundaryType == NEWMAN)) {
      datamid[0] = data[0];
      return datamid[0];
    }
    else if ((j==0.5) && (boundaryType== DIRICHLET))
      return datamid[0];

    if ((j==n_points+0.5) && (boundaryType == NEWMAN)) {
      datamid[n_points] = data[n_points-1];
      return datamid[n_points];
    }
    else if ((j==n_points+0.5) && (boundaryType == DIRICHLET))
      return datamid[n_points];
    
    i = (int)floorf(j);
    if ( isFlux == 1 )
      return datamid[i];

    datamid[i] = 0.5*( (*this)[i] + (*this)[i+1] );
    return datamid[i];
  }
  else if ( (checkInt != 0.5) && (checkInt != 0) ){
    cout << "ERROR: fractional index is not plus minus 0.5"<<endl;
    cout << "index = "<<checkInt<<endl;
  }
  else {
    cout << "ERROR: index out range"<<endl;
    exit(1);
  }
}


double &gridFunc1D::operator[]( float j )
{
  int i;
  float checkInt = j - floorf(j);

  // check if j is integer
  if (checkInt == 0) {
    i=(int)j;
    if ( (1<=i) && (i<=n_points) )
      return data[i-1];
    // for i higher than N
    else if ( (n_points<i) && (i<2*n_points) && (boundaryType == 0) )
      return data[i-n_points-1];
    // for i lower than 1
    else if ( (-n_points<i) && (i<1) && (boundaryType == 0) )
      return data[n_points-i-1];
    else {
      cout << "ERROR: index out range"<<endl;
      cout << "Index requested = "<<i<<endl;
      exit(1);
    }
  }

  // check if j is not integer
  else if ( checkInt == 0.5 ){
    // check if we are at the boundaries
    if ((j==0.5) && (boundaryType == NEWMAN)) {
      datamid[0] = data[0];
      return datamid[0];
    }
    else if ((j==0.5) && (boundaryType== DIRICHLET))
      return datamid[0];

    if ((j==n_points+0.5) && (boundaryType == NEWMAN)) {
      datamid[n_points] = data[n_points-1];
      return datamid[n_points];
    }
    else if ((j==n_points+0.5) && (boundaryType == DIRICHLET))
      return datamid[n_points];
    
    i = (int)floorf(j);
    if ( isFlux == 1 )
      return datamid[i];

    datamid[i] = 0.5*( (*this)[i] + (*this)[i+1] );
    return datamid[i];
  }
  else if ( (checkInt != 0.5) && (checkInt != 0) ){
    cout << "ERROR: fractional index is not plus minus 0.5"<<endl;
    cout << "index = "<<checkInt<<endl;
  }
  else {
    cout << "ERROR: index out range"<<endl;
    exit(1);
  }
}


//! Get the number of points
int gridFunc1D::Npoints() const
{
  return n_points;
}


gridFunc1D gridFunc1D::operator+( const gridFunc1D &B ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] + B.data[i];

  return C;
}


gridFunc1D gridFunc1D::operator-( const gridFunc1D &B ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] - B.data[i];

  return C;
}

gridFunc1D operator-( const double a, const gridFunc1D &B )
{
  gridFunc1D C( B.n_points );

  for( int i=0; i<B.n_points; i++ )
    C.data[i] = a - B.data[i];

  return C;
}


gridFunc1D gridFunc1D::operator*( const gridFunc1D &B ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] * B.data[i];

  return C;
}


gridFunc1D gridFunc1D::operator*( const double &b ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] * b;

  return C;
}


gridFunc1D operator*( const double &a, const gridFunc1D &B )
{
  gridFunc1D C( B.n_points );

  for( int i=0; i<B.n_points; i++ )
    C.data[i] = a * B.data[i];

  return C;
}


gridFunc1D gridFunc1D::operator/( const double &b ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] / b;

  return C;
}


gridFunc1D gridFunc1D::operator/( const gridFunc1D &B ) const
{
  gridFunc1D C( n_points );

  for( int i=0; i<n_points; i++ )
    C.data[i] = data[i] / B.data[i];

  return C;
}


gridFunc1D operator/( const double a, const gridFunc1D &B )
{
  gridFunc1D C( B.n_points );

  for( int i=0; i<B.n_points; i++ )
    C.data[i] = a / B.data[i];

  return C;
}


const gridFunc1D &gridFunc1D::operator=( const gridFunc1D &B )
{
  // data may be NULL when using the mVector class template. A vector
  // of gridFunc1D may be created, but the gridFunc1D constructor is never
  // caller with the appropiate size parameter. To avoid segmentation
  // errors we check here to give it the righ space.
  if ( data == NULL )
    create( B.n_points );

  if ( n_points != B.n_points ){
    cerr << endl;
    cerr << "ERROR: operator= gridFunc1D's are not the same size" << endl;
    exit(1);
  }
  

  for( int i=0; i<n_points; i++ )
    data[i] = B.data[i];

  return *this;
}


const gridFunc1D &gridFunc1D::operator=( const double b )
{
  if ( data == NULL ){
    cerr << endl;
    cerr << "ERROR: operator= gridFunc1D does't have allocated memory" << endl;
    exit(1);
  }

  for( int i=0; i<n_points; i++ )
    data[i] = b;

  return *this;
}

const gridFunc1D &gridFunc1D::operator=( const int b )
{
  if ( data == NULL ){
    cerr << endl;
    cerr << "ERROR: operator= gridFunc1D does't have allocated memory" << endl;
    exit(1);
  }

  for( int i=0; i<n_points; i++ )
    data[i] = (double)b;

  return *this;
}

gridFunc1D sqrt( const gridFunc1D &A )
{
  gridFunc1D temp( A.Npoints() );

  for( int i=0; i<A.Npoints(); i++ )
    temp.data[i] = sqrt( A.data[i] );

  return temp;
}

//! Checks is a given value is greater than a
void gridFunc1D::ishuge( double a )
{
  for( int i=1; i<=n_points; i++ ){
    if ( fabs(data[i]) > a )
      cerr << "Point "<<i<<" = "<<data[i]<<endl;
  }
}


ostream &operator<<( ostream &os, const gridFunc1D &A )
{
  os << "gridFunc1D:" << endl;

  for( int i=1; i<=A.Npoints(); i++ )
    os << "["<<i<<"] = " << A[i] << endl;

  return os;
}


