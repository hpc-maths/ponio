#include <iostream>
#include <valarray>
#include <numeric>
#include <numbers>

#include <solver/problem.hpp>
#include <solver/solver.hpp>
#include <solver/observer.hpp>
#include <solver/butcher_methods.hpp>

/*
solve Arenstorf orbit problem
*/
  
struct arenstorf_model
{
    using state_t = std::valarray<double>;
  
    double mu;
  
    arenstorf_model(double m) : mu(m) {}
   
    state_t operator () (double t , state_t const& y)
    { 
        double y1 = y[0], y2 = y[1], y3 = y[2], y4 = y[3];
        double r1 = sqrt((y1+mu)*(y1+mu) + y2*y2);
        double r2 = sqrt((y1-1+mu)*(y1-1+mu) + y2*y2);
        double dy1 = y3;
        double dy2 = y4;
        double dy3 = y1 + 2*y4 - (1-mu)*(y1+mu)/(r1*r1*r1) - mu*(y1-1+ mu)/(r2*r2*r2);
        double dy4 = y2 - 2*y3 - (1-mu)*y2/(r1*r1*r1) - mu*y2/(r2*r2*r2);
        return {dy1, dy2, dy3, dy4};
    }
};

int main(int, char**)
{
    std::string dirname = "arenstorf_data";  
    std::filesystem::create_directories(dirname);
    std::string filename;

    using state_t = std::valarray<double>;

    double tf = 17.0652165601579625588917206249;
    double dt = 1e-5;

    double mu = 0.012277471;

    auto arenstorf_pb = ode::make_problem(arenstorf_model(mu));

    state_t yini = { 0.994, 0., 0., -2.00158510637908252240537862224 };

    filename = std::filesystem::path(dirname) / "arenstorf_rk546m.dat";
    ode::solve(arenstorf_pb, ode::butcher::rk54_6m<>(1e-5), yini, {0.,tf}, dt, observer::file_observer(filename));

    filename = std::filesystem::path(dirname) / "arenstorf_rk547m.dat";
    ode::solve(arenstorf_pb, ode::butcher::rk54_7m<>(1e-5), yini, {0.,tf}, dt, observer::file_observer(filename));

    filename = std::filesystem::path(dirname) / "arenstorf_rk547s.dat";
    ode::solve(arenstorf_pb, ode::butcher::rk54_7s<>(1e-5), yini, {0.,tf}, dt, observer::file_observer(filename));

    return 0;
}