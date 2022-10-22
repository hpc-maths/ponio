#include <iostream>
#include <valarray>
#include <sstream>
#include <vector>
#include <string>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

/*
Lotka-Volterra system
---------------------

$$
  \begin{cases}
    \frac{\mathrm{d}x}{\mathrm{d}t} = \alpha x - \beta xy\\
    \frac{\mathrm{d}y}{\mathrm{d}t} = \delta xy - \gamma y\\
  \end{cases}
$$

with $\alpha = \frac{2}{3}$, $\beta = \frac{4}{3}$, and $\delta = \gamma = 1$. Initial condition $(x,y)(t=0) = (x_0,x_0)$ done by user.

This system is solved by RK(11,8) Runge-Kutta method with time step $\Delta t=0.1$ to the time $t\leq 15$.
 */

int main(int argc, char** argv)
{
     // default filename
    std::string dirname = "lv_data";  
    std::filesystem::create_directories(dirname);
    auto filename = std::filesystem::path(dirname) / "lv.dat";

    using state_t = std::valarray<double>;

    double x0 = 1.0; // default x0
    if (argc > 1) {
        filename = argv[1];
        x0       = std::stof(argv[2]);
    }
    observer::file_observer fobs(filename);

    double alpha=2./3., beta=4./3., gamma=1., delta=1.; // parameter
    auto lotka_volterra_pb = ode::make_simple_problem( // define problem
            [=]( double t , state_t const& u ) -> state_t {
                return {
                    alpha*u[0] - beta*u[0]*u[1] ,
                    delta*u[0]*u[1] - gamma*u[1]
                };
            }
        );

    std::vector<double> t_span = {0.,15.}; // begin and end time
    double dt = 0.1; // time step
    state_t u0 = {x0, x0}; // initial condition
    ode::solve(lotka_volterra_pb, ode::butcher::rk_118<>(), u0, t_span, dt, fobs);

    return 0;
}