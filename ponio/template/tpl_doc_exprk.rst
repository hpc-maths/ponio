List of exponential Runge-Kutta methods
=======================================

{% for rk in list_exprk %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endfor %}
