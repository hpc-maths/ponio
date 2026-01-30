List of explicit Runge-Kutta methods
====================================

{% for rk in list_erk %}{% if rk.b2 is undefined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}
