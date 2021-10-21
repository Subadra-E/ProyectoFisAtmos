#include <iomanip>
#include "gridFunc1D.hpp"
#include "Euler1D.hpp"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int sgn( double val )
{
  return (0. < val) - (val < 0.);
}


//! Namespace to hold all necesary constants
namespace Eu1D
{
  const double gam = 5.0/3.0;   //! gamma, ratio between Cp and Cv 
}

//!
//! q is a mVector that hold the variables
//! q[0] is the density
//! q[1] is velocity
//! q[2] is pressure 
//!
Euler1D::Euler1D( const parameterReader &p ) : q(3), r_m(3), r_0(3), r_p(3), alpha(3), lambda(3), fluxRHS(3), Flux(3)
{
  const int n  = atoi( p.getParam( "n_points" ) );
  const int cf = atoi( p.getParam( "convergence_factor" ) );
 
  // allocate memory
  moment.create( n );
  energ.create( n );
  
  for( int j=0; j<3; j++ ){
    q[j].create( n );
    Flux[j].create( n );
    lambda[j].create( n );
    r_m[j].create( n );
    r_0[j].create( n );
    r_p[j].create( n );
    alpha[j].create( n );
    fluxRHS[j].create( n );
  }
  
  convFactor = cf;
}

mVector<gridFunc1D> &Euler1D::stateVector()
{
  return q;
}

const gridFunc1D &Euler1D::density() const
{
  return q[0];
}

const gridFunc1D &Euler1D::vel() const
{
  return q[1];
}

const gridFunc1D &Euler1D::pres() const
{
  return q[2];
}

const gridFunc1D &Euler1D::momentum() 
{
  moment = q[0]*q[1];
  return moment;
}

const gridFunc1D &Euler1D::energy()
{
  energ = 0.5*q[0]*q[1]*q[1] + q[2]/(Eu1D::gam-1.);
  return energ;
}



// Pointwise formulas
//! Momentum in terms of density and velocity
inline double momentum_PW( double d, double v )
{
  return d*v;
}

inline double energy_PW( double d, double v, double p )
{
  return 0.5*d*v*v + p/(Eu1D::gam-1);
}

//! Four types of shock tubes as initial data
void Euler1D::initialData( gridFunc1D &x, const parameterReader &p )
{
  const int N = q[0].Npoints();
  const std::string s( p.getParam( "initial_data_type" ) );

  if ( s == "shock_tube_1" ){
    for( int i=1; i<=N; i++ ){
      // shock tube 1
      q[1][i] = 0.0;
      q[2][i] = (i<N/2)? 3.0 : 1.0;
      q[0][i] = (i<N/2)? 3.0 : 1.0;
    }
  } else if ( s == "shock_tube_2" ){
    for( int i=1; i<=N; i++ ){
      // shock tube 2
      q[1][i] = 0.0;
      q[2][i] = (i<N/2)? 10.0 : 1.0;
      q[0][i] = 1.0;
    }
  } else if ( s == "shock_tube_3" ){
    for( int i=1; i<=N; i++ ){
      // shock tube 3
      q[1][i] = (i<N/2)? 3.0 : 1.0;
      q[2][i] = 1.0;
      q[0][i] = (i<N/2)? 1 : 2.0;
    }
  } else if ( s == "shock_tube_4" ){
    for( int i=1; i<=N; i++ ){
      // shock tube 4
      q[1][i] = (i<N/4)? 0.0 : 1.0;
      q[2][i] = 1.0;
      q[0][i] = (i<N/4)? 1.0 : 1.0;
    }
  } else if ( s == "gaussian1" ){
    for( int i=1; i<=N; i++ ){
      // gaussian
      q[0][i] = 3.0 + exp( -50*pow(x[i]-x[N]/2, 2) );
      q[1][i] = 0.0;
      q[2][i] = 3.0;
    }
  } else if ( s == "gaussian2" ){
    for( int i=1; i<=N; i++ ){
      // gaussian
      q[0][i] = 3.0;
      q[1][i] = exp( -50*pow(x[i]-x[N]/2, 2) );
      q[2][i] = 3.0;
    }
  }  else if ( s == "gaussian3" ){
    for( int i=1; i<=N; i++ ){
      // gaussian
      q[0][i] = 3.0;
      q[1][i] = 0.0;
      q[2][i] = 3.0 + exp( -50*pow(x[i]-x[N]/2, 2) );
    }
  } else if ( s == "whatever" ){
    for( int i=1; i<=N; i++ ){

      // densidad
      q[0][i] = exp(x[i]);

      // velocidad
      q[1][i] = 0.0;

      // presion
      q[2][i] = 4.0 + sin(M_PI*x[i]);
    }
  }
  else {
    std::cerr << endl;
    std::cerr << "ERROR: Value "<<s<<" of parameter initial_data_type is invalid"<<endl;
    std::cerr << "Valid options are: shock_tube_1"<<endl
              << "                   shock_tube_2"<<endl
              << "                   shock_tube_3"<<endl
              << "                   shock_tube_4"<<endl
              << "                   gaussian1"<<endl
	      << "                   gaussian2"<<endl
	      << "                   gaussian3"<<endl
	      << "                   whatever"<<endl;
    exit( 1 );
  }  
}

//! Evolve in time from t to t+dt
//!
//! This is the most simple thing we can do. It is
//! a first order method
//!
mVector<gridFunc1D> Euler1D::RHS( const mVector<gridFunc1D> &iq, double dx )
{
  // iq stands for input q
  const int N = iq[0].Npoints();
  
  calcFlux3( iq, alpha, r_m, r_0, r_p, lambda, Flux );
  
  // compute net fluxes
  for( int i=3; i<=N-2; i++ ){
    fluxRHS[0][i] = -( Flux[0][i+1] - Flux[0][i] )/dx;
    fluxRHS[1][i] = -( Flux[1][i+1] - Flux[1][i] )/dx;
    fluxRHS[2][i] = -( Flux[2][i+1] - Flux[2][i] )/dx;
  }

  return fluxRHS;
}


void boundaryCondition( mVector<gridFunc1D> &q )
{
  const int N = q[0].Npoints();
  // constant extrapolation at boundaries
  q[0][1] = q[0][2] = q[0][3];
  q[1][1] = q[1][2] = q[1][3];
  q[2][1] = q[2][2] = q[2][3];
  
  q[0][N] = q[0][N-1] = q[0][N-2];
  q[1][N] = q[1][N-1] = q[1][N-2];
  q[2][N] = q[2][N-1] = q[2][N-2];
}






void Euler1D::output( gridFunc1D &x, double t )
{
  static int count = 0;
  ofstream out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, out11, out12, out13, out14, out15, out16, out17, out18, out19, out20;
  std::ios_base::openmode mode;

  if ( count == 0 )
    mode = ios::out;
  else
    mode = ios::app;

  out1.open( "density.dat",  mode );
  out2.open( "momentum.dat", mode );
  out3.open( "energy.dat",   mode ); 
  out4.open( "vel.dat",      mode );
  out5.open( "pres.dat",     mode );
  /* 
     out6.open( "densityConv.dat",     mode );
  out7.open( "flux0rhs.dat", mode );
  out8.open( "flux1rhs.dat", mode );
  out9.open( "flux2rhs.dat", mode );
  out10.open( "flux0.dat", mode );
  out11.open( "flux1.dat", mode );
  out12.open( "flux2.dat", mode );
  out13.open( "velConv.dat", mode );
  out14.open( "presConv.dat", mode );
  out15.open( "flux0rhsConv.dat", mode );
  out16.open( "flux1rhsConv.dat", mode );
  out17.open( "flux2rhsConv.dat", mode );
  out18.open( "flux0Conv.dat", mode );
  out19.open( "flux1Conv.dat", mode );
  out20.open( "flux2Conv.dat", mode );
  */

  density() .outputGnuplot( out1, x, t );
  momentum().outputGnuplot( out2, x, t );
  energy()  .outputGnuplot( out3, x, t );
  vel()     .outputGnuplot( out4, x, t );
  pres()    .outputGnuplot( out5, x, t );
  // fluxRHS[0].outputByColumn( out7, x, t );
  // fluxRHS[1].outputByColumn( out8, x, t );
  // fluxRHS[2].outputByColumn( out9, x, t );
  // Flux[0].outputByColumn( out10, x, t );
  // Flux[1].outputByColumn( out11, x, t );
  // Flux[2].outputByColumn( out12, x, t );

  // output for convergence test
  /*
  gridFunc1D xConv( q[0].Npoints()/convFactor );
  gridFunc1D densityConv( q[0].Npoints()/convFactor );
  gridFunc1D velConv( q[0].Npoints()/convFactor );
  gridFunc1D presConv( q[0].Npoints()/convFactor );
  gridFunc1D flux0rhsConv( q[0].Npoints()/convFactor );
  gridFunc1D flux1rhsConv( q[0].Npoints()/convFactor );
  gridFunc1D flux2rhsConv( q[0].Npoints()/convFactor );
  gridFunc1D flux0Conv( q[0].Npoints()/convFactor );
  gridFunc1D flux1Conv( q[0].Npoints()/convFactor );
  gridFunc1D flux2Conv( q[0].Npoints()/convFactor );
  convergenceOutput( x, q[0], xConv, densityConv, convFactor );
  convergenceOutput( x, q[1], xConv, velConv, convFactor );
  convergenceOutput( x, q[2], xConv, presConv, convFactor );
  convergenceOutput( x, fluxRHS[0], xConv, flux0rhsConv, convFactor );
  convergenceOutput( x, fluxRHS[1], xConv, flux1rhsConv, convFactor );
  convergenceOutput( x, fluxRHS[2], xConv, flux2rhsConv, convFactor );
  convergenceOutput( x, Flux[0], xConv, flux0Conv, convFactor );
  convergenceOutput( x, Flux[1], xConv, flux1Conv, convFactor );
  convergenceOutput( x, Flux[2], xConv, flux2Conv, convFactor );
  

  densityConv.outputByColumn( out6, xConv, t );
  velConv.outputByColumn( out13, xConv, t );
  presConv.outputByColumn( out14, xConv, t );
  flux0rhsConv.outputByColumn( out15, xConv, t );
  flux1rhsConv.outputByColumn( out16, xConv, t );
  flux2rhsConv.outputByColumn( out17, xConv, t );
  flux0Conv.outputByColumn( out18, xConv, t );
  flux1Conv.outputByColumn( out19, xConv, t );
  flux2Conv.outputByColumn( out20, xConv, t );
  */


  out1.close();
  out2.close();
  out3.close();
  out4.close();
  out5.close();
  /*
  out7.close();
  out8.close();
  out9.close();
  out10.close();
  out11.close();
  out12.close();
  out13.close();
  out14.close();
  out15.close();
  out17.close();
  out18.close();
  out19.close();
  out20.close();
  */

  
  count++;
}


void convergenceOutput( gridFunc1D &x, gridFunc1D &dens, gridFunc1D &xc, gridFunc1D &densc, int convFactor )
{
  int n = x.Npoints() / convFactor;

  if ( convFactor == 2 ){
    for( int i=1; i<=n; i++ ){
      //xc[i] = 0.5*(x[2*i-1] + x[2*i] );
      //densc[i] = 0.5*(dens[2*i-1] + dens[2*i] );
      xc[i] = x[2*i-1];
      densc[i] = dens[2*i-1];
    }
  }
  else if ( convFactor == 4 ){
    for( int i=1; i<=n; i++ ){
      //xc[i] = 0.25*(x[4*i-3] + x[4*i-2] + x[4*i-1] + x[4*i] );
      //densc[i] = (-dens[4*i-3] + 9*dens[4*i-2] + 9*dens[4*i-1] - dens[4*i] )/16.0;
      xc[i] = x[4*i-3];
      densc[i] = dens[4*i-3];
    }
  }
}




//!
//! Compute the flux f(q) of the equation system based on the cell
//! averages
//!
//void Euler1D::calc_fluxRHS()
//{
//  fluxRHS[0] = q[1]*q[0];
//  fluxRHS[1] = q[1]*q[1]*q[0] + q[2];
//  fluxRHS[2] = q[1] * ( energy() + q[2] );
//}






double minmod( double a, double b )
{
  if ( a*b < 0.0 ) return 0.0;
  else if ( fabs(a) < fabs(b) ) return a;
  else return b;
}

double superbee( double theta )
{
  return fmax(fmax(0, fmin(1, 2*theta)), fmin(2, theta));
}


void calc_dq_version2( mVector<gridFunc1D> &qL, mVector<gridFunc1D> &qR, 
	      mVector<gridFunc1D> &alpha, 
	      mVector<gridFunc1D> &r_m, 
	      mVector<gridFunc1D> &r_0, 
	      mVector<gridFunc1D> &r_p,
	      mVector<gridFunc1D> &lambda )
{
  double u, c, pL, pR, c2, u2, h, hL, hR, srhoL, srhoR, rho, EL, ER;
  double dq[3], v0[3], v1[3], v2[3];
   
  using Eu1D::gam;

  for( int i=3; i<=qL[0].Npoints()-1; i++ ){
    srhoL = sqrt(qL[0][i]);
    srhoR = sqrt(qR[0][i]);
    rho = sqrt( qL[0][i]*qR[0][i] );
    //rho = 0.5*(qL[0][i] + qR[0][i] );
    u  = (srhoL*qL[1][i] + srhoR*qR[1][i])/(srhoL + srhoR);
    //u = 0.5*(qL[1][i] + qR[1][i] );
    u2 = u*u;
    EL = 0.5*qL[0][i]*qL[1][i]*qL[1][i] + qL[2][i]/(gam-1.0);
    ER = 0.5*qR[0][i]*qR[1][i]*qR[1][i] + qR[2][i]/(gam-1.0);
    hL = ( EL + qL[2][i] ) / qL[0][i];
    hR = ( ER + qR[2][i] ) / qR[0][i];
    h  = (srhoL*hL + srhoR*hR)/(srhoL + srhoR);
    //h = 0.5*(hL + hR);
    c2 = (gam-1.0)*(h - 0.5*u2 );
    c  = sqrt(c2);

    dq[0] = qR[0][i] - qL[0][i];
    dq[1] = qR[1][i] - qL[1][i];
    dq[2] = qR[2][i] - qL[2][i];

    // matrix components to get alpha in conserved variables
    /*
    v0[0] = 0.25*( 2.*u/c +    u2/c2*(gam-1.) );
    v0[1] = 0.25*( 4.     - 2.*u2/c2*(gam-1.) );
    v0[2] = 0.25*(-2.*u/c +    u2/c2*(gam-1.) );
    
    v1[0] = 0.25*(-2./c - 2.*u/c2*(gam-1.) );
    v1[1] = 0.25*(        4.*u/c2*(gam-1.) );
    v1[2] = 0.25*( 2./c - 2.*u/c2*(gam-1.) );
    
    v2[0] = 0.25*( 2./c2*(gam-1.) );
    v2[1] = 0.25*(-4./c2*(gam-1.) );
    v2[2] = 0.25*( 2./c2*(gam-1.) );
    */
    // matrix components to get alpha in primitive variables
    v0[0] = 0.;
    v0[1] = 1.;
    v0[2] = 0.;

    v1[0] =-0.5*rho/c;
    v1[1] = 0.;
    v1[2] = 0.5*rho/c;
  
    v2[0] = 0.5/c2;
    v2[1] =-1.0/c2; 
    v2[2] = 0.5/c2;
    
  
    for( int j=0; j<3; j++ )
      alpha[j][i] = v0[j]*dq[0] + v1[j]*dq[1] + v2[j]*dq[2];

    // eigenvector for system in conserved variables
    /*
    r_m[0][i] = 1;
    r_m[1][i] = u - c;
    r_m[2][i] = h - u*c;

    r_0[0][i] = 1;
    r_0[1][i] = u;
    r_0[2][i] = 0.5*u2;

    r_p[0][i] = 1;
    r_p[1][i] = u + c;
    r_p[2][i] = h + u*c;
    */

    // eigenvector for system in primitive variables
    
    r_m[0][i] = 1;
    r_m[1][i] = -c/rho;
    r_m[2][i] = c*c;

    r_0[0][i] = 1.;
    r_0[1][i] = 0.;
    r_0[2][i] = 0.;

    r_p[0][i] = 1;
    r_p[1][i] = c/rho;
    r_p[2][i] = c*c;
    


    lambda[0][i] = u - c;
    lambda[1][i] = u;
    lambda[2][i] = u + c;
  }
}


double calc_qL( double q_im2, double q_im1, double q_i )
{
  double dq1, dq2, dq;

  dq1 = q_i - q_im1;
  dq2 = q_im1 - q_im2;
  dq = minmod( dq1, dq2 );

  return q_im1 + 0.5*dq;
}

double calc_qR( double q_im1, double q_i, double q_ip1 )
{
  double dq1, dq2, dq;

  dq1 = q_ip1 - q_i;
  dq2 = q_i - q_im1;
  dq = minmod( dq1, dq2 );

  return q_i - 0.5*dq;
}


void calc_flux( mVector<gridFunc1D> &q, mVector<gridFunc1D> &f )		  
{
  const int N = q[0].Npoints();
  gridFunc1D energy(N);

  energy = 0.5*q[0]*q[1]*q[1] + q[2]/(Eu1D::gam-1.0);
  f[0] = q[1]*q[0];
  f[1] = q[1]*q[1]*q[0] + q[2];
  f[2] = q[1] * ( energy + q[2] );

}


void calcFlux3( const mVector<gridFunc1D> &iq, 
		mVector<gridFunc1D> &alpha,
		mVector<gridFunc1D> &r_m,
		mVector<gridFunc1D> &r_0,
		mVector<gridFunc1D> &r_p,
		mVector<gridFunc1D> &lambda,
		mVector<gridFunc1D> &Flux )
{
  mVector<gridFunc1D> qR(3), qL(3), fluxL(3), fluxR(3);
  const int N = iq[0].Npoints();

  for( int j=0; j<3; j++ ){
    qL[j].create( N );
    qR[j].create( N );
    fluxL[j].create( N );
    fluxR[j].create( N );
  }

  for( int i=3; i<=N-1; i++ ){
    for( int j=0; j<3; j++ ){
      qL[j][i] = calc_qL( iq[j][i-2], iq[j][i-1], iq[j][i]   );
      qR[j][i] = calc_qR( iq[j][i-1], iq[j][i],   iq[j][i+1] );
    }
  }

  calc_dq_version2( qL, qR, alpha, r_m, r_0, r_p, lambda );

  calc_flux( qL, fluxL );
  calc_flux( qR, fluxR );

  for( int i=3; i<=N-1; i++ ){
    for( int j=0; j<3; j++ ){
      Flux[j][i] = 
	0.5*(fluxL[j][i]+fluxR[j][i]) 
	- 0.5*( alpha[0][i]*fabs(lambda[0][i])*r_m[j][i] +
		alpha[1][i]*fabs(lambda[1][i])*r_0[j][i] +
		alpha[2][i]*fabs(lambda[2][i])*r_p[j][i] );
      
      
    

    }
    /*
    cout << "L["<<i<<"] = "<<setprecision(10)<<fluxL[0][i]<<endl;
    cout << "R["<<i<<"] = "<<setprecision(10)<<fluxR[0][i]<<endl;
    cout << "alpha[]  = "<< alpha[0][i]<<" "<<alpha[1][i]<<" "<<alpha[2][i]<<endl;
    cout << "lambda[] = "<< lambda[0][i]<<" "<<lambda[1][i]<<" "<<lambda[2][i]<<endl;
    cout << "r_m[]    = "<< r_m[0][i]<<" "<<r_m[1][i]<<" "<<r_m[2][i]<<endl;
    cout << "r_0[]    = "<< r_0[0][i]<<" "<<r_0[1][i]<<" "<<r_0[2][i]<<endl;
    cout << "r_p[]    = "<< r_p[0][i]<<" "<<r_p[1][i]<<" "<<r_p[2][i]<<endl;
    */

  }
}
 









//========================================================================================


/*
void Euler1D::calc_q_eigen()
{
  const int n = density().Npoints();
  gridFunc1D u(n), c(n), c2(n), u2(n), h(n);
  double q[3];
  mVector<gridFunc1D> v0(3), v1(3), v2(3);
  // memory for v0, v1, v2 is allocated by operator=

  using Eu1D::gam;
 
  //calc_vel();
  //calc_pres();

  // this gives the fluxes f(q)
  calc_fluxRHS();

  u = vel;
  u2 = vel*vel;
  c2 = gam*pres/density();
  c  = sqrt(c2);
  h = (energy() + pres)/density();

  v0[0] = 0.25*( 2.*u/c +    u2/c2*(gam-1.) );
  v0[1] = 0.25*( 4.     - 2.*u2/c2*(gam-1.) );
  v0[2] = 0.25*(-2.*u/c +    u2/c2*(gam-1.) );

  v1[0] = 0.25*(-2./c - 2.*u/c2*(gam-1.) );
  v1[1] = 0.25*(        4.*u/c2*(gam-1.) );
  v1[2] = 0.25*( 2./c - 2.*u/c2*(gam-1.) );
  
  v2[0] = 0.25*( 2./c2*(gam-1.) );
  v2[1] = 0.25*(-4./c2*(gam-1.) );
  v2[2] = 0.25*( 2./c2*(gam-1.) );

  // components of q
  alpha = q[0]*v0 + q[1]*v1 + q[2]*v2;

  // component of f(q)
  phi = fluxRHS[0]*v0 + fluxRHS[1]*v1 + fluxRHS[2]*v2;

  // It is important to put the decimal dot. Otherwize it is interpreted as int
  r_m[0] = 1.;
  r_m[1] = u - c;
  r_m[2] = h - u*c;

  r_0[0] = 1.;
  r_0[1] = u;
  r_0[2] = 0.5*u2;

  r_p[0] = 1.;
  r_p[1] = u + c;
  r_p[2] = h + u*c;

  lambda[0] = u - c;
  lambda[1] = u;
  lambda[2] = u + c;

  //cout << lambda[0] <<" "<< lambda[1] <<" "<< lambda[2] << endl;
 
}
 */

/*
double calc_phi_m( double l_i, double l_im1, double phi, double L, double omega )
{
  if ( (l_i > 0.) && (l_im1 > 0.) )
    return 0.;
  else if ( (l_i < 0. ) && (l_im1 < 0. ) )
    return phi;
  else
    return 0.5*( phi - L*omega );
}
*/

/*
double calc_phi_p( double l_i, double l_im1, double phi, double L, double omega )
{
  if ( (l_i < 0.) && (l_im1 < 0.) )
    return 0.;
  else if ( (l_i > 0. ) && (l_im1 > 0. ) )
    return phi;
  else
    return 0.5*( phi - L*omega );
}
*/

/*
void Euler1D::MarquinaFlux()
{
  double Lambda[3], phi_m[3], phi_p[3];

  for( int i=2; i<=q[0].Npoints(); i++ ){
    Lambda[0] = fmax( lambda[0][i-1], lambda[0][i] );
    Lambda[1] = fmax( lambda[1][i-1], lambda[1][i] );
    Lambda[2] = fmax( lambda[2][i-1], lambda[2][i] );

    for( int j=0; j<3; j++ ){
      phi_m[j] = calc_phi_m( lambda[j][i], lambda[j][i-1], phi[j][i],   Lambda[j], alpha[j][i] );
      phi_p[j] = calc_phi_p( lambda[j][i], lambda[j][i-1], phi[j][i-1], Lambda[j], alpha[j][i-1] );
    }

    for( int j=0; j<3; j++ ){
      MarqFlux[j][i] = 
	phi_p[0]*r_m[j][i-1] + phi_m[0]*r_m[j][i] + 
	phi_p[1]*r_0[j][i-1] + phi_m[1]*r_0[j][i] +
	phi_p[2]*r_p[j][i-1] + phi_m[2]*r_p[j][i];

      Flux_m[j][i] = MarqFlux[j][i] - fluxRHS[j][i-1];
      Flux_p[j][i] = fluxRHS[j][i] - MarqFlux[j][i];
    }
  }
}
*/


// flux and Roe solver to 1st order
/*
void Euler1D::calcFlux2( double dt, double dx )
{
  calc_fluxRHS();

  for( int i=2; i<=density().Npoints(); i++ ){
    calc_dq( i, q, alpha, r_m, r_0, r_p, lambda );

    Flux[0][i] = 
      0.5*(fluxRHS[0][i-1]+fluxRHS[0][i]) 
      - 0.5*( alpha[0][i]*fabs(lambda[0][i])*r_m[0][i] +
	      alpha[1][i]*fabs(lambda[1][i])*r_0[0][i] +
	      alpha[2][i]*fabs(lambda[2][i])*r_p[0][i] );

    Flux[1][i] = 
      0.5*(fluxRHS[1][i-1]+fluxRHS[1][i]) 
      - 0.5*( alpha[0][i]*fabs(lambda[0][i])*r_m[1][i] +
	      alpha[1][i]*fabs(lambda[1][i])*r_0[1][i] +
	      alpha[2][i]*fabs(lambda[2][i])*r_p[1][i] );

    Flux[2][i] = 
      0.5*(fluxRHS[2][i-1]+fluxRHS[2][i]) 
      - 0.5*( alpha[0][i]*fabs(lambda[0][i])*r_m[2][i] +
	      alpha[1][i]*fabs(lambda[1][i])*r_0[2][i] +
	      alpha[2][i]*fabs(lambda[2][i])*r_p[2][i] );
  }
}  
*/

/*
void calc_dq( int i, const mVector<gridFunc1D> &q, 
	      mVector<gridFunc1D> &alpha, 
	      mVector<gridFunc1D> &r_m, 
	      mVector<gridFunc1D> &r_0, 
	      mVector<gridFunc1D> &r_p,
	      mVector<gridFunc1D> &lambda )
{
  double u, c, pL, pR, c2, u2, h, hL, hR, srhoL, srhoR, rho;
  double dq[3], v0[3], v1[3], v2[3];
  
  //u  = mom[i]/den[i];
  //u2 = u*u;
  //p  = (gam-1.0)*(ene[i] - 0.5*den[i]*u2 );
  //c2 = gam*p/den[i];
  //c  = sqrt(c2);
  //h = (ene[i] + p)/den[i];
  
  
  using Eu1D::gam;

  srhoL = sqrt(q[0][i-1]);
  srhoR = sqrt(q[0][i]);
  //rho = sqrt(den[i-1]*den[i]);

  u  = (srhoL*q[1][i-1]/q[0][i-1] + srhoR*q[1][i]/q[0][i])/(srhoL + srhoR);
  u2 = u*u;
  pL = (gam-1.0)*(q[2][i-1] - 0.5*q[0][i-1]*pow(q[1][i-1]/q[0][i-1],2) );
  pR = (gam-1.0)*(q[2][i  ] - 0.5*q[0][i  ]*pow(q[1][i  ]/q[0][i  ],2) );
  hL = (q[2][i-1] + pL ) / q[0][i-1];
  hR = (q[2][i  ] + pR ) / q[0][i  ];
  h  = (srhoL*hL + srhoR*hR)/(srhoL + srhoR);
  c2 = (gam-1.0)*(h - 0.5*u2 );
  c  = sqrt(c2);


  dq[0] = q[0][i] - q[0][i-1];
  dq[1] = q[1][i] - q[1][i-1];
  dq[2] = q[2][i] - q[2][i-1];
  
  v0[0] = 0.25*( 2*u/c +   u2/c2*(gam-1) );
  v0[1] = 0.25*( 4     - 2*u2/c2*(gam-1) );
  v0[2] = 0.25*(-2*u/c +   u2/c2*(gam-1) );

  v1[0] = 0.25*(-2/c - 2*u/c2*(gam-1) );
  v1[1] = 0.25*(       4*u/c2*(gam-1) );
  v1[2] = 0.25*( 2/c - 2*u/c2*(gam-1) );
  
  v2[0] = 0.25*( 2/c2*(gam-1) );
  v2[1] = 0.25*(-4/c2*(gam-1) );
  v2[2] = 0.25*( 2/c2*(gam-1) );
  
  for( int j=0; j<3; j++ )
    alpha[j][i] = v0[j]*dq[0] + v1[j]*dq[1] + v2[j]*dq[2];

  r_m[0][i] = 1;
  r_m[1][i] = u - c;
  r_m[2][i] = h - u*c;

  r_0[0][i] = 1;
  r_0[1][i] = u;
  r_0[2][i] = 0.5*u2;

  r_p[0][i] = 1;
  r_p[1][i] = u + c;
  r_p[2][i] = h + u*c;

  lambda[0][i] = u - c;
  lambda[1][i] = u;
  lambda[2][i] = u + c;
}
*/


/*
void Euler1D::calcFlux( double dt, double dx )
{
  mVector<gridFunc1D> lam_p(3), lam_m(3);
  
  lam_p[0].create( q[0].Npoints() );
  lam_p[1].create( q[0].Npoints() );
  lam_p[2].create( q[0].Npoints() );
  lam_m[0].create( q[0].Npoints() );
  lam_m[1].create( q[0].Npoints() );
  lam_m[2].create( q[0].Npoints() );


  for( int i=2; i<=q[0].Npoints(); i++ ){
    calc_dq( i, q, alpha, r_m, r_0, r_p, lambda );

    for( int j=0; j<3; j++ ){
      lam_p[j][i] = lam_m[j][i] = 0.0;
      lambda[j][i] > 0.0 ? lam_p[j][i] = lambda[j][i] : lam_m[j][i] = lambda[j][i];
    }
  }
     
  // Aplus
  Flux_p = lam_p[0]*alpha[0]*r_m + lam_p[1]*alpha[1]*r_0 + lam_p[2]*alpha[2]*r_p;
  // Aminus
  Flux_m = lam_m[0]*alpha[0]*r_m + lam_m[1]*alpha[1]*r_0 + lam_m[2]*alpha[2]*r_p;
}
*/

/*
void Euler1D::calcLWFlux( double dt, double dx )
{
  int I;
  double W_i_m[3], W_i_0[3], W_i_p[3];
  double W_I_m[3], W_I_0[3], W_I_p[3];
  double mag_i_m, mag_i_0, mag_i_p, mag_I_m, mag_I_0, mag_I_p;
  double theta_m, theta_0, theta_p, phi_m, phi_0, phi_p;
  double lam0, lam1, lam2;
  double F_LW_i[3];

  for( int i=3; i<=q[0].Npoints()-1; i++ ){
    W_i_m[0] = alpha[0][i]*r_m[0][i] - alpha[0][i-1]*r_m[0][i-1];
    W_i_m[1] = alpha[0][i]*r_m[1][i] - alpha[0][i-1]*r_m[1][i-1];
    W_i_m[2] = alpha[0][i]*r_m[2][i] - alpha[0][i-1]*r_m[2][i-1];

    W_i_0[0] = alpha[1][i]*r_0[0][i] - alpha[1][i-1]*r_0[0][i-1];
    W_i_0[1] = alpha[1][i]*r_0[1][i] - alpha[1][i-1]*r_0[1][i-1];
    W_i_0[2] = alpha[1][i]*r_0[2][i] - alpha[1][i-1]*r_0[2][i-1];

    W_i_p[0] = alpha[2][i]*r_p[0][i] - alpha[2][i-1]*r_p[0][i-1];
    W_i_p[1] = alpha[2][i]*r_p[1][i] - alpha[2][i-1]*r_p[1][i-1];
    W_i_p[2] = alpha[2][i]*r_p[2][i] - alpha[2][i-1]*r_p[2][i-1];

    lam0 = 0.5*( lambda[0][i-1] + lambda[0][i] );
    lam0 > 0.0 ? I=i-1 : I=i+1;
    W_I_m[0] = alpha[0][I]*r_m[0][I] - alpha[0][I-1]*r_m[0][I-1];
    W_I_m[1] = alpha[0][I]*r_m[1][I] - alpha[0][I-1]*r_m[1][I-1];
    W_I_m[2] = alpha[0][I]*r_m[2][I] - alpha[0][I-1]*r_m[2][I-1];

    lam1 = 0.5*( lambda[1][i-1] + lambda[1][i] );
    lam1 > 0.0 ? I=i-1 : I=i+1;
    W_I_0[0] = alpha[1][I]*r_0[0][I] - alpha[1][I-1]*r_0[0][I-1];
    W_I_0[1] = alpha[1][I]*r_0[1][I] - alpha[1][I-1]*r_0[1][I-1];
    W_I_0[2] = alpha[1][I]*r_0[2][I] - alpha[1][I-1]*r_0[2][I-1];

    lam2 = 0.5*( lambda[2][i-1] + lambda[2][i] );
    lam2 > 0.0 ? I=i-1 : I=i+1;
    W_I_p[0] = alpha[2][I]*r_p[0][I] - alpha[2][I-1]*r_p[0][I-1];
    W_I_p[1] = alpha[2][I]*r_p[1][I] - alpha[2][I-1]*r_p[1][I-1];
    W_I_p[2] = alpha[2][I]*r_p[2][I] - alpha[2][I-1]*r_p[2][I-1];

    // calc magnitude for unitary vector in the direction of W_i
    mag_i_m = sqrt( pow(W_i_m[0],2) + pow(W_i_m[1],2) + pow(W_i_m[2],2) );
    mag_i_0 = sqrt( pow(W_i_0[0],2) + pow(W_i_0[1],2) + pow(W_i_0[2],2) );
    mag_i_p = sqrt( pow(W_i_p[0],2) + pow(W_i_p[1],2) + pow(W_i_p[2],2) );


    // project W_I onto W_i
    //mag_I_m = (W_I_m[0]*W_i_m[0] + W_I_m[1]*W_i_m[1] + W_I_m[2]*W_i_m[2])/mag_i_m;
    //mag_I_0 = (W_I_0[0]*W_i_0[0] + W_I_0[1]*W_i_0[1] + W_I_0[2]*W_i_0[2])/mag_i_0;
    //mag_I_p = (W_I_p[0]*W_i_p[0] + W_I_p[1]*W_i_p[1] + W_I_p[2]*W_i_p[2])/mag_i_p;

    //theta_m = mag_I_m/mag_i_m;
    //theta_0 = mag_I_0/mag_i_0;
    //theta_p = mag_I_p/mag_i_p;
    
    phi_m = 1;//minmod( 1.0, theta_m );
    phi_0 = 1;//minmod( 1.0, theta_0 );
    phi_p = 1;//minmod( 1.0, theta_p );
   
    
    //phi_m = superbee( theta_m );
    //phi_0 = superbee( theta_0 );
    //phi_p = superbee( theta_p );
    
    //cout << mag_i_m <<" "<< mag_i_0 <<" "<< mag_i_p << endl;
    //cout << W_i_p[0] <<" "<< W_i_p[1] <<" "<< W_i_p[2] << endl;

    // calc Lax-Wendroff flux term
    for( int j=0; j<3; j++ )
      F_LW_i[j] = 
	0.5*( fabs(lam0) * ( 1.0 - dt/dx*fabs(lam0) ) * phi_m*W_i_m[j] +
	      fabs(lam1) * ( 1.0 - dt/dx*fabs(lam1) ) * phi_0*W_i_0[j] +
	      fabs(lam2) * ( 1.0 - dt/dx*fabs(lam2) ) * phi_p*W_i_p[j] );

    Flux_LW[0][i] = F_LW_i[0];
    Flux_LW[1][i] = F_LW_i[1];
    Flux_LW[2][i] = F_LW_i[2];
    
    
    //if ( i==100 )
    //  cout << "phi_m   = " << phi_m 
	   << " phi_0  = " << phi_0 
	   << " phi_p  = " << phi_p  << endl
	   << " mag_i_m = " << mag_i_m
	   << " mag_i_0 = " << mag_i_0
	   << " mag_i_p = " << mag_i_p << endl;
    

  }
  
  //cout << Flux_LW[0] <<" "<< Flux_LW[1] <<" "<< Flux_LW[2] << endl;

}
*/
