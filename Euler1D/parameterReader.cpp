
#include "parameterReader.hpp"

//! Constructor
//!
//! In the constructor the Parameters file is read. The pairs "parameter name"
//! and "parameter value" are stored in a matrix.
//!
parameterReader::parameterReader( std::string fileName )
{
  std::ifstream infile( fileName.c_str() );
  std::string line, paramName, paramValue;
  std::size_t pos;
  int c = 0;

  while( std::getline( infile, line ) ){
    // get position of the = sign
    pos = line.find( "=" );
    // get the string before and after the = sign
    paramName  = line.substr( 0, pos-1 );
    paramValue = line.substr( pos+1 );
    // remove space infront of parameter value
    paramValue.erase(remove_if(paramValue.begin(), paramValue.end(), isspace),paramValue.end());
    // assing to paramTable
    paramTable[c][0] = paramName;
    paramTable[c][1] = paramValue;

    //cout << paramTable[c][0] << "qqq" << paramTable[c][1] << "q"<<endl;
    c++;
  }

  len = c;
  infile.close();
}

//! Returns the value of the parameter name given.
//!
//! The function takes the parameter name s and returns a the value as a const char*
//! If the value is numeric, it has to be converted to int or double.
//!
const char *parameterReader::getParam( std::string s ) const
{
  int i=0;

  while( paramTable[i][0] != s ){
    i++;
    if ( i>len ){
      std::cerr << "ERROR: Parameter "<<s<<" was not found in the parameter file"<<endl;
      exit(1);
    }
  }

  return paramTable[i][1].c_str();
}



/* test main

int main()
{

  parameterReader p( "Parameters" );

  cout << p.getParam( "n_points" ) << endl;

  return 0;
}

*/
