List of embedded Runge-Kutta methods
====================================

Explicit methods
~~~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}

Diagonal implicit methods
~~~~~~~~~~~~~~~~~~~~~~~~~

{% for rk in list_dirk %}{% if rk.b2 is defined %}
.. doxygentypedef:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}

Lawson methods
~~~~~~~~~~~~~~

{% for rk in list_erk %}{% if rk.b2 is defined %}
.. doxygentypedef:: ponio::runge_kutta::l{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}
