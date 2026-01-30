List of Lawson Runge-Kutta methods
==================================

{% for rk in list_erk %}{% if rk.b2 is undefined %}
.. doxygenvariable:: ponio::runge_kutta::l{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygenvariable:: ponio::runge_kutta::l{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}
