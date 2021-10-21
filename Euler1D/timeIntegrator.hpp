#ifndef timeIntegrator_hpp
#define timeIntegrator_hpp

#include "Euler1D.hpp"

class TimeIntegrator
{
public:
  TimeIntegrator();
void integrate( Euler1D &, double, double );
};

#endif
