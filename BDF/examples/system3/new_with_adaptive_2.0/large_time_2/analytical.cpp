#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>

using namespace std;
using namespace boost::numeric::odeint;

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
ofstream myfile;
void write_solution( const vector_type &x , const double t )
{  
    
    myfile << t << '\t' << x[0] << '\t' << x[1] <<'\t' << x[2] << endl;
    cout << t << '\t' << x[0] << '\t' << x[1] <<'\t' << x[2] <<endl;
    
}


int main(int argc, char **argv)
{
      
      vector_type x(3);
      x(0) = 1, x(1) = 0, x(2) = 0;

      myfile.open( "analytical.txt");
      integrate( stiff_system() , x , 0.0 , 2700.0 , 100.0 , write_solution );  
      myfile.close(); 
}
