#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "bdf2.hpp"
#include <fstream>

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{    
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ]; 
        dxdt[ 1 ] = x[ 0 ];
    }
};

struct stiff_system_jacobi
{
     void operator()( const vector_type & /* x */ , matrix_type &J , const double & /* t */ )
    {
        J( 0 , 0 ) = -101.0;
        J( 0 , 1 ) = -100.0;
        J( 1 , 0 ) = 1.0;
        J( 1 , 1 ) = 0.0;
    }
};

using namespace std;

int main( int argc , char** argv )
{  
    
    ofstream myfile;
    myfile.open( "bdf2.txt");

    vector_type x(2);
    x(0) = 2, x(1) = 1;
    double t = 0.0;
    double dt = 1E-4;// initial conditions
    int steps = 0;
    
    boost::numeric::odeint::bdf2< vector_type , double > solver;
   
    while ( steps < 15000 )
      {  
          solver.do_step( make_pair( stiff_system() , stiff_system_jacobi() ), x , t , dt );
          t += dt;
        { 
           cout<<"  "<< t <<"  " << x[0] << "  "<< x[1] <<"\n";
           myfile << t << "  " << x[0] << "  " << x[1] << "\n";
        }
          steps++;
      }
    myfile.close();
}
