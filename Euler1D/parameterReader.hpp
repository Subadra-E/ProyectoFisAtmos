#ifndef parameterReader_H
#define parameterReader_H

#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

using std::endl;
using std::cout;

class parameterReader
{
private:
  // This would read a parameter file of a thousand lines. Hopefully we will not
  // need more.
  int len; //!< the length of the parameter table 
  std::string paramTable[1000][2]; //!< a len by 2 matrix holding the parameter name and the parameter value

public:
  parameterReader( std::string );
  const char *getParam( std::string s ) const;
};

#endif



