#include <iostream>
#include <tuple>
#include <filesystem>

#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

// solve $\dot{u} = u$ with $u(t=0) = 1$, and $t\in[0,2]$.

int main (int, char**)
{
    std::string dirname = "exp_data";  
    std::filesystem::create_directories(dirname);
    auto filename = std::filesystem::path(dirname) / "exp.dat";
    observer::file_observer fobs(filename);

    auto identity = [](double t, double u){ return u; };
    double x0 = 1.0;
    double dt = 0.1;

    ode::solve(identity, ode::butcher::rk_nssp_21<>(), x0, {0.,2.0}, dt, fobs);

    return 0;
    }