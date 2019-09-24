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
        dxdt[ 0 ] = -0.04 * x[ 0 ] + 10000.0 * x[ 1 ] * x[ 2 ];
        dxdt[ 1 ] =  0.04 * x[ 0 ] - 10000.0 * x[ 1 ] * x[ 2 ] - 3 * 10000000.0 * x[ 1 ] * x[ 1 ];
        dxdt[ 2 ] =  3 * 10000000.0 * x[ 1 ] * x[ 1 ];
    }   
};

struct stiff_system_jacobi
{
     void operator()( const vector_type & x , matrix_type &J , const double & /* t */ )
    {
        J( 0 , 0 ) = -0.04;
        J( 0 , 1 ) = 10000.0 * x [2];
        J( 0 , 2 ) = 10000.0 * x[1];
        J( 1 , 0 ) = 0.04;
        J( 1 , 1 ) = -10000.0 * x[2] - 3 * 10000000.0 * 2 * x[1];
        J( 1 , 2 ) = -10000.0 * x[1];
        J( 2 , 0 ) = 0;
        J( 2 , 1 ) = 3 * 10000000.0 * 2 * x[1];
        J( 2 , 2 ) = 0;
    }
};

using namespace std;

int main( int argc , char** argv )
{  
    
    vector_type x(3);
    x(0) = 1.0, x(1) = 0.0, x(2) = 0.0;
    double t = 0.0;
    double dt = 0.0001;//initial values

    ofstream myfile;
    myfile.open( "bdf2.txt");
    boost::numeric::odeint::bdf2< boost::numeric::ublas::vector<double>, double > solver;
   
    for( size_t i=0 ; i< 100 ; ++i , t+= dt )
     {
            solver.do_step( make_pair( stiff_system() , stiff_system_jacobi() ), x , t , dt ); 
            cout << t << " " << x[0] << " " << x[1] <<"  "<< x[2] <<"\n";    
            myfile << t << " " << x[0] << " " << x[1] <<"  "<< x[2] <<"\n";     
     }
    myfile.close();
    return 0;
}
