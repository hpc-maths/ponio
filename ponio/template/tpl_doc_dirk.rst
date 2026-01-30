List of diagonal implicit Runge-Kutta methods
=============================================

{% for rk in list_dirk %}{% if rk.b2 is undefined %}
.. doxygenfunction:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

{% for rk in list_dirk %}{% if rk.b2 is defined %}
.. doxygenfunction:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}
