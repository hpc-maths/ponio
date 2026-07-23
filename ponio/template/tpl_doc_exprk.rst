List of exponential Runge-Kutta methods
=======================================

{% for rk in list_exprk %}

{{ rk.label }}
~~~~~~~~~~~~~~

+ **name:** {{ rk.label }}
+ **label in ponio:** :cpp:type:`ponio::runge_kutta::{{ rk.id }}_t`
+ **stages:** {{ rk.A|length }}
+ **order:** {{ rk.order }}
{% if 'bib' in rk %}+ **bibliography:** [{{ rk.bib.bib }}]({{ rk.bib.url }}){%- endif %}

.. math::

    \begin{array}{c|{%- for ci in rk.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.butcher.A %}
          {{ rk.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.butcher.b|join(' & ') }}
    \end{array}

{% endfor %}
