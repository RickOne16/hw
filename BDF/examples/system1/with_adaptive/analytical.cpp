#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double > vector_type;

struct stiff_system
{    
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ]; 
        dxdt[ 1 ] = x[ 0 ];
    }
};
ofstream myfile;
void write_solution( const vector_type &x , const double t )
{  
    
    myfile << t << '\t' << x[0] << '\t' << x[1] << endl;
    cout << t << '\t' << x[0] << '\t' << x[1] << endl;
    
}


int main(int argc, char **argv)
{
      
      vector_type x(2);
      x(0) = 2E-4, x(1) = 1E-4;

      myfile.open( "analytical.txt");
      integrate( stiff_system() , x , 0.0 , 50.0 , 0.01 , write_solution );  
      myfile.close(); 
}
