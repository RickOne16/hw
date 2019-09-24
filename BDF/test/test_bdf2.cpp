#include <boost/config.hpp>
#ifdef BOOST_MSVC
    #pragma warning(disable:4996)
#endif

#define BOOST_TEST_MODULE odeint_bdf2

#include <utility>
#include <iostream>

#include "bdf2.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::unit_test;
using namespace boost::numeric::odeint;

typedef double value_type;
typedef boost::numeric::ublas::vector< value_type > vector_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;

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

BOOST_AUTO_TEST_SUITE( bdf2_test )

BOOST_AUTO_TEST_CASE( test_bdf2_stepper )

{
    typedef bdf2< state_type, value_type > stepper_type;
    stepper_type stepper;
           
    state_type x(2);
    x(0) = 2.0, x(1) = 1.0;
   
    double t = 0;
    double dt = 0.0001;

    stepper.do_step( std::make_pair(stiff_system(), stiff_system_jacobi() ), x, t, dt);


}


BOOST_AUTO_TEST_SUITE_END()


