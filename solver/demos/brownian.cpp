#include <iostream>
#include <valarray>
#include <numeric>
#include <random>
#include <sstream>
#include <filesystem>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

int main (int argc, char** argv)
{
    std::string dirname = "brownian_data";

    std::size_t n = 10;
    if (argc > 1) { n = std::atoi(argv[1]); }

    using state_t = std::valarray<double>;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<> d{0.,2};

    double dt = 1e-3;

    auto brownian_pb = ode::make_simple_problem(
        [&](double t, state_t const& y)->state_t {
            return {d(gen), d(gen)};
        }
    );

    state_t yini = { 0., 0. };

    for ( auto i=0u ; i<n ; ++i ) 
    {
        std::stringstream ssfilename; ssfilename << "brownian_"<< i << ".dat";
        auto filename = std::filesystem::path(dirname) / ssfilename.str();
        observer::file_observer fobs(filename);
        ode::solve(brownian_pb, ode::butcher::rk_33(), yini, {0.,10.}, dt, fobs);
    }

  return 0;
}