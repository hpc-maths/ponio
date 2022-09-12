#include <doctest/doctest.h>

#include <solver/butcher_methods.hpp>
#include <solver/solver.hpp>

#include "compute_order.hpp"

{% for rk in list_meth %} 
TEST_CASE("order::{{ rk.id }}")
{
  using rk_t = ode::butcher::{{ rk.id }}<>;
  CHECK( order<rk_t>() == doctest::Approx(rk_t::order).epsilon(0.05) );
}
{% endfor %}

