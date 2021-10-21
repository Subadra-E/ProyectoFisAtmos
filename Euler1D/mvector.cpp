
#include "mvector.hpp"
#include "gridFunc1D.hpp"

//template void mVector<double>::set( int i, double temp );
//template mVector<double> operator*( double a, const mVector<double> &z );

template mVector<gridFunc1D> operator*( double a, const mVector<gridFunc1D> &z );
template mVector<gridFunc1D> operator*( const gridFunc1D &a, const mVector<gridFunc1D> &z );


//! Constructor
template <class T>
mVector<T>::mVector()
{
  dim = 0;
  vec = NULL;
}


template <class T>
mVector<T>::mVector( int d )
{
  vec = NULL;
  resize(d);
}

template< class T >
mVector<T>::~mVector()
{
  erase();
}


// There is no initialization because it is a template
template <class T>
void mVector<T>::resize( int d )
{
  erase();

  dim = d;
  vec = new T[d];
}


template <class T>
void mVector<T>::erase()
{
  if ( vec != NULL )
    delete[] vec;

  vec = NULL;
}

//! Copy constructor
template<class T>
mVector<T>::mVector( const mVector<T> &z )
{
  vec = NULL;
  resize( z.dim );

  // now we copy also the components
  for( int i=0; i<dim; i++ )
    vec[i] = z.vec[i];  
}


//! Set the value of each components
template <class T>
void mVector<T>::set( int i, T temp )
{
  if ( (0<=i) && (i<dim) )
    vec[i] = temp;
  else{
    cerr << "vector::set index out of range" << endl;
    exit(1);
  }
}


template<class T>
int mVector<T>::getDim() const
{
  return dim;
}


//! Overload operator[] to get the values of each component
template<class T>
T &mVector<T>::operator[]( int i ) const
{
  return vec[i];
}


//! Overload operator[] to set the values of each component
template<class T>
T &mVector<T>::operator[]( int i )
{
  return vec[i];
}



//----------------------------------
//----- overloading operators ------
//----------------------------------

template<class T>
const mVector<T> &mVector<T>::operator=( const mVector<T> &z )
{
  dim = z.dim;

  for( int i=0; i<dim; i++ )
    vec[i] = z.vec[i];  

  return *this;
}


template<class T>
const mVector<T> &mVector<T>::operator=( const double a )
{
  if ( vec == NULL ){
    cerr << endl;
    cerr << "ERROR: mvector objecto does not have allocated memory" << endl;
    exit(1);
  }

  for( int i=0; i<dim; i++ )
    vec[i] = a;  

  return *this;
}

template<class T>
const mVector<T> &mVector<T>::operator=( const int a )
{
  if ( vec == NULL ){
    cerr << endl;
    cerr << "ERROR: mvector objecto does not have allocated memory" << endl;
    exit(1);
  }

  for( int i=0; i<dim; i++ )
    vec[i] = (double)a;  

  return *this;
}


template<class T>
mVector<T> mVector<T>::operator+( const mVector<T> &z ) const
{
  // temp will have non-initialize elements. The operator= for the
  // class of element must assign the proper information.
  mVector<T> temp(dim);

  for( int i=0; i<dim; i++ )
    temp.vec[i] = vec[i] + z.vec[i];

  return temp;
}


template<class T>
mVector<T> mVector<T>::operator-( const mVector<T> &z ) const
{
  mVector<T> temp(dim);

  for( int i=0; i<dim; i++ )
    temp.vec[i] = vec[i] - z.vec[i];

  return temp;
}


template<class T>
mVector<T> mVector<T>::operator*( const double a ) const
{
  mVector<T> temp(dim);

  for( int i=0; i<dim; i++ )
    temp.vec[i] = vec[i] * a;

  return temp;
}


template<class T>
mVector<T> operator*( const double a, const mVector<T> &z )
{
  mVector<T> temp( z.getDim() );

  for( int i=0; i<z.getDim(); i++ )
    temp[i] = a * z[i];

  return temp;
}

template<class T>
mVector<T> operator*( const gridFunc1D &a, const mVector<T> &z )
{
  mVector<T> temp( z.getDim() );

  for( int i=0; i<z.getDim(); i++ )
    temp[i] = a * z[i];

  return temp;
}






//! L2 norm of the vector
/*
template <class T>
double mVector<T>::L2norm() const
{
  double sum = 0.0;

  for( int i=0; i<getDim(); i++ )
    sum += pow( vec[i], 2 );

  return sqrt( sum );
}
*/	 
	 


template class mVector<gridFunc1D>;

/* test code

int main()
{
  mVector<double> A(4), B(4);

  A.set( 0, 10.0 );
  A.set( 1, 20.0 );
  A.set( 2, 30.0 );
  A.set( 3, 35.0 );

  B = A*4;

  cout << A.L2norm() << "  " << B.L2norm() << endl;


  A[0] = 3;
  A[1] = 4;
  A[2] = 0;
  A[3] = 0;

  cout << A[0] << " " << A[1] << " " << A[2] << endl;
  cout << A.L2norm() << endl;




  return 0;
}


*/

