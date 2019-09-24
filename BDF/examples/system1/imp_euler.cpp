#include <boost/numeric/odeint/stepper/implicit_euler.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct stiff_system
{    
    void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
    {
        dxdt[ 0 ] = -101.0 * x[ 0 ] - 100.0 * x[ 1 ]; //Had there been a coefficient of 't' then dfdt != 0
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
    clock_t begin, end;
    double time_spent;
    vector_type x(2);
    x(0) = 2.0, x(1) = 1.0;
    double t = 0.0;
    double dt = 0.0001;
    begin = clock();
    boost::numeric::odeint::implicit_euler<double> solver;
    for( size_t i=0 ; i< 301 ; ++i , t+= dt )
    {
        solver.do_step( make_pair( stiff_system() , stiff_system_jacobi() ), x , t , dt );
        cout << t << " " << x[0] << " " << x[1] << "\n";
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    cout << " time spent : " << time_spent<<"seconds"<< endl;
    return 0;
}
