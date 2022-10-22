#include <iostream>
#include <valarray>
#include <numeric>
#include <vector>
#include <filesystem>

#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

int main(int, char**)
{
    std::string dirname = "lorenz_data";  
    std::filesystem::create_directories(dirname);
    auto filename = std::filesystem::path(dirname) / "lorenz.dat";
    observer::file_observer fobs(filename);

    using state_t = std::valarray<double>;

    double sigma=10. , rho=28., beta=8./3.;
    auto lorenz = [=](double t, state_t const& u) -> state_t
    {
        auto du1 =  sigma*(u[1] - u[0]);
        auto du2 =  rho*u[0] - u[1] - u[0]*u[2];
        auto du3 = u[0]*u[1] - beta*u[2];
        return {du1, du2, du3};
    };
  
    state_t u0 = {1.,1.,1.};
    std::vector<double> tspan = {0.,20.};
    double dt = 0.01;

    ode::solve(lorenz, ode::butcher::rk_nssp_53<>(), u0, tspan, dt, fobs);

    return 0;
}