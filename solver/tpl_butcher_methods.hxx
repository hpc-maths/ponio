#pragma once

#include "butcher_tableau.hpp"
#include "generic_butcher_rk.hpp"

namespace ode::butcher {

{% for rk in list_meth %}
template <typename value_t=double>
struct butcher_{{ rk.label }} : public butcher_tableau<{{ rk.A|length }},value_t>
{
  using base_t = butcher_tableau<{{ rk.A|length }},value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_{{ rk.label }}()
  : base_t(
    {{ '{{' }}
        {% for ai in rk.A -%}
        { {{ ai }} }{{ "," if not loop.last else "" }}
        {%- endfor %}
    {{ '}}' }}, // A
    { {{ rk.b }} }, // b
    { {{ rk.c }} }  // c
  )
  {}
};
template <typename value_t=double>
using {{ rk.label }} = runge_kutta::explicit_rk_butcher<butcher_{{ rk.label }}<value_t>>;

{% endfor %}

} // namespace ode::butcher
