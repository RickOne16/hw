#ifndef BOOST_NUMERIC_ODEINT_STEPPER_BDF2_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_STEPPER_BDF2_HPP_INCLUDED

#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/stepper/implicit_euler.hpp>

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/operations_dispatcher.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint/stepper/detail/rotating_buffer.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace boost {
namespace numeric {
namespace odeint {




template< class state , class ValueType , class Algebra = typename algebra_dispatcher< state >::algebra_type, class Operations = typename operations_dispatcher< state >::operations_type ,class Resizer = initially_resizer >
class bdf2
{


public:

    typedef ValueType time_type;
    typedef ValueType value_type;
    typedef state state_type ;
    typedef state_wrapper< state_type > wrapped_state_type;
    typedef state_type deriv_type;
    typedef state_wrapper< deriv_type > wrapped_deriv_type;
    typedef boost::numeric::ublas::matrix< ValueType > matrix_type;
    typedef state_wrapper< matrix_type > wrapped_matrix_type;
    typedef boost::numeric::ublas::permutation_matrix< size_t > pmatrix_type;
    typedef state_wrapper< pmatrix_type > wrapped_pmatrix_type;
    typedef Resizer resizer_type;
    typedef stepper_tag stepper_category;
    typedef bdf2< state , ValueType , Algebra, Operations, Resizer > stepper_type;

    typedef detail::rotating_buffer< state_type, 2> rotbuffer_type;
    rotbuffer_type m_previous_step;    

   // bool flag = true;

public:

    bdf2( value_type epsilon = 1E-6 , bool flag = true): m_epsilon( epsilon ), m_flag( flag )
    { }

    template< typename System >
    void do_step( System system , state_type& x , time_type t , time_type dt )
    {
        
        
        typedef typename odeint::unwrap_reference< System >::type system_type;
        typedef typename odeint::unwrap_reference< typename system_type::first_type >::type deriv_func_type;
        typedef typename odeint::unwrap_reference< typename system_type::second_type >::type jacobi_func_type;
        system_type &sys = system;
        deriv_func_type &deriv_func = sys.first;
        jacobi_func_type &jacobi_func = sys.second;   
        
        m_resizer.adjust_size( x , detail::bind( &stepper_type::template resize_impl<state_type> , detail::ref( *this ) , detail::_1 ) );
        
            
        for( size_t i=0 ; i<x.size() ; ++i )
            m_pm.m_v[i] = i;

        if( m_flag )
        {  
            m_previous_step[0] = x;

            implicit_euler< double > imp_euler;
            imp_euler.do_step(system, x, t, dt);

            m_previous_step[1] = x;
            t += dt;//t = t(k+1)
            m_flag = false;
           
        }
          
        else
        {
        
            t += dt; // t = t(k+2)

            deriv_func( m_previous_step[1] , m_dxdt.m_v , t );// calculating x'(k+2)
            m_b.m_v = (2.0/3.0) * dt * m_dxdt.m_v;
        
            jacobi_func( m_previous_step[1] , m_jacobi.m_v  , t );
            m_jacobi.m_v *= dt;
            m_jacobi.m_v -= boost::numeric::ublas::identity_matrix< value_type >( x.size() );
        
            solve( m_b.m_v , m_jacobi.m_v );
            m_x.m_v = (4.0/3.0) * m_previous_step[1]- (1.0/3.0) * m_previous_step[0]- m_b.m_v; // BDF2 formula
        
       
            while( boost::numeric::ublas::norm_2( m_b.m_v ) > m_epsilon )//Newton's iteration for accuracy
            {
                deriv_func( m_x.m_v , m_dxdt.m_v , t );
                m_b.m_v = (4.0/3.0) * m_previous_step[1] - (1.0/3.0) * m_previous_step[0] - m_x.m_v + (2.0/3.0) * dt * m_dxdt.m_v;
                solve( m_b.m_v , m_jacobi.m_v );
                m_x.m_v -= m_b.m_v;     
            }

            m_previous_step.rotate();
            m_previous_step[1] = m_x.m_v; 
            x = m_x.m_v;

        }

    }// do_step
        
private:

    template< class StateIn >
    bool resize_impl( const StateIn &x )
    {
        bool resized = false;
        resized |= adjust_size_by_resizeability( m_dxdt , x , typename is_resizeable<deriv_type>::type() );
        resized |= adjust_size_by_resizeability( m_x , x , typename is_resizeable<state_type>::type() );
        resized |= adjust_size_by_resizeability( m_b , x , typename is_resizeable<deriv_type>::type() );
        resized |= adjust_size_by_resizeability( m_jacobi , x , typename is_resizeable<matrix_type>::type() );
        resized |= adjust_size_by_resizeability( m_pm , x , typename is_resizeable<pmatrix_type>::type() );
        return resized;
    }


    void solve( state_type &x , matrix_type &m )
    {
        int res = boost::numeric::ublas::lu_factorize( m , m_pm.m_v );
        if( res != 0 ) exit(0);
        boost::numeric::ublas::lu_substitute( m , m_pm.m_v , x );
    }

private:

    value_type m_epsilon;
    resizer_type m_resizer;
    wrapped_deriv_type m_dxdt;
    wrapped_state_type m_x;
    wrapped_deriv_type m_b;
    wrapped_matrix_type m_jacobi;
    wrapped_pmatrix_type m_pm;
    bool m_flag;

};
} // odeint
} // numeric
} // boost

#endif

