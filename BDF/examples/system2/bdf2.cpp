#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "bdf2.hpp"


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{    
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -80.6 * x[ 0 ] + 119.4 * x[ 1 ]; //Had there been a coefficient of 't' then dfdt != 0
        dxdt[ 1 ] = 79.6 * x[ 0 ] - 120.4 * x[ 1 ];
    }
};

struct stiff_system_jacobi
{
     void operator()( const vector_type & /* x */ , matrix_type &J , const double & /* t */ )
    {
        J( 0 , 0 ) = -80.6;
        J( 0 , 1 ) = 119.4;
        J( 1 , 0 ) = 79.6;
        J( 1 , 1 ) = -120.4;
    }
};

using namespace std;

int main( int argc , char** argv )
{  
    
    vector_type x(2);
    x(0) = 2.0, x(1) = 1.0;//initial values
    double t = 0.0;
    double dt = 0.0001;
    
    boost::numeric::odeint::bdf2< boost::numeric::ublas::vector<double>, double > solver;
   
    for( size_t i=0 ; i< 301; ++i , t+= dt )
     {
         solver.do_step( make_pair( stiff_system() , stiff_system_jacobi() ), x , t , dt );
          cout << t << " " << x[0] << " " << x[1] << "\n";          
     }
    
    return 0;
}
