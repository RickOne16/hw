#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <time.h>

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
     void operator()( const vector_type & /* x */ ,matrix_type &J , const double & /* t */ ,vector_type &dfdt)
    {
        J( 0 , 0 ) = -101.0;
        J( 0 , 1 ) = -100.0;
        J( 1 , 0 ) = 1.0;
        J( 1 , 1 ) = 0.0;
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
    }
};
using namespace std;

int main( int argc , char** argv )
{  
    clock_t begin, end;
    double time_spent;

    vector_type x(2);
    vector_type xout(2), xerr(2);

    xout(0) = 0.0, xout(1)=0.0;
    x(0) = 2.0, x(1) = 1.0;
    double t = 0.0;
    double dt = 0.0001;//initial conditions

    for ( int j=0; j< 500; j++)
  {

    time_spent = 0.0;

    vector_type x(2);
    x(0) = 2.0, x(1) = 1.0;

    double t = 0.0;
    double dt = 0.0001;//initial conditions

    begin = clock();    

    boost::numeric::odeint::rosenbrock4<double> solver;
    for( size_t i=0 ; i< 301 ; ++i , t+= dt )
    {
        solver.do_step( make_pair( stiff_system() , stiff_system_jacobi() ), x , t , xout ,dt, xerr );
        x = xout;
    }

    end = clock();

    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    cout << j << "  " << time_spent<< endl;
  }
    return 0;
}
