#include <cmath>
#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/time_span.hpp>

/*
solve Curtiss and Hirschfelder problem:

$$
    \begin{aligned}
        \dot{y} =  k(\cos(t) - y) \\
        y(0) = y_0
    \end{aligned}
$$

$y_0 = 2$
*/

int
main()
{
    using namespace ponio::observer;

    double const k   = 50.;
    double const y_0 = 2.;

    ponio::time_span<double> t_span = { 0., 4.0 };
    double dt                       = 0.05;

    { // explicit Runge-Kutta method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::rk_33(), y_0, t_span, dt, "ch_erk.txt"_fobs );
    }

    { // embedded Runge-Kutta method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::rk54_6m().abs_tol( 1e-6 ).rel_tol( 1e-4 ), y_0, t_span, dt, "ch_dp.txt"_fobs );
    }

    { // diagonal implicit Runge-Kutta method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };
        auto df = [=]( double t, double /* y */ )
        {
            return -k;
        };

        auto pb = ponio::make_implicit_problem( f, df );

        ponio::solve( pb, ponio::runge_kutta::dirk34(), y_0, t_span, dt, "ch_dirk.txt"_fobs );
    }

    { // Lawson method
        double L = -k;
        auto N   = [=]( double t, double y, double& dy )
        {
            dy = k * std::cos( t );
        };

        auto pb = ponio::make_lawson_problem( L, N );

        auto my_exp = []( double x )
        {
            return std::exp( x );
        };

        ponio::solve( pb, ponio::runge_kutta::lrk_33( my_exp ), y_0, t_span, dt, "ch_lrk.txt"_fobs );
    }

    { // exponential Runge-Kutta method
        double L = -k;
        auto N   = [=]( double t, double y, double& dy )
        {
            dy = k * std::cos( t );
        };

        auto pb = ponio::make_lawson_problem( L, N );

        ponio::solve( pb, ponio::runge_kutta::exprk22(), y_0, t_span, dt, "ch_exprk.txt"_fobs );
    }

    { // Runge-Kutta Chebyshev method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::explicit_rkc2<5>(), y_0, t_span, dt, "ch_rkc.txt"_fobs );
    }

    { // ROCK2 method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::rock::rock2(), y_0, t_span, dt, "ch_rock2.txt"_fobs );
    }

    { // ROCK4 method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::rock::rock2(), y_0, t_span, dt, "ch_rock4.txt"_fobs );
    }

    { // Runge-Kutta Legendre first-order method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::explicit_rkl1<5>(), y_0, t_span, dt, "ch_rkl1.txt"_fobs );
    }

    { // Runge-Kutta Legendre second-order method
        auto f = [=]( double t, double y, double& dy )
        {
            dy = k * ( std::cos( t ) - y );
        };

        ponio::solve( f, ponio::runge_kutta::explicit_rkl2<5>(), y_0, t_span, dt, "ch_rkl2.txt"_fobs );
    }

    { // Lie splitting method
        auto f1 = [=]( double t, double y, double& dy )
        {
            dy = k * std::cos( t );
        };
        auto f2 = [=]( double t, double y, double& dy )
        {
            dy = -k * y;
        };

        auto pb = ponio::make_problem( f1, f2 );

        auto lie = ponio::splitting::make_lie_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ),
            std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ) );

        ponio::solve( pb, lie, y_0, t_span, dt, "ch_lie.txt"_fobs );
    }

    { // Strang splitting method
        auto f1 = [=]( double t, double y, double& dy )
        {
            dy = k * std::cos( t );
        };
        auto f2 = [=]( double t, double y, double& dy )
        {
            dy = -k * y;
        };

        auto pb = ponio::make_problem( f1, f2 );

        auto strang = ponio::splitting::make_strang_tuple( std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ),
            std::make_pair( ponio::runge_kutta::rk_33(), 0.5 * dt ) );

        ponio::solve( pb, strang, y_0, t_span, dt, "ch_strang.txt"_fobs );
    }
}
