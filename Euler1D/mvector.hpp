#ifndef mvector_hpp
#define mvectop_hpp

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "gridFunc1D.hpp"


using namespace std;


template <class T>
class mVector
{
private:
  int dim;
  T *vec;

public:
  mVector();
  ~mVector();
  mVector( int );
  void resize( int );
  mVector( const mVector<T> & );
  void erase();
  void set( int, T );
  int getDim() const;
  T &operator[]( int i ) const;
  T &operator[]( int i );
  const mVector<T> &operator=( const mVector<T> & );
  const mVector<T> &operator=( const double a );
  const mVector<T> &operator=( const int a );
  mVector<T> operator+( const mVector<T> & ) const;
  mVector<T> operator-( const mVector<T> & ) const;
  mVector<T> operator*( const double a ) const;
  //double L2norm() const;

};

template <class T>
mVector<T> operator*( const double a, const mVector<T> &z );
template <class T>
mVector<T> operator*( const gridFunc1D &a, const mVector<T> &z );





#endif
