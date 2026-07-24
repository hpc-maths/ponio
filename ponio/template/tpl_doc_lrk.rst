List of Lawson Runge-Kutta methods
==================================

{% for rk in list_erk %}{% if rk.b2 is undefined %}

Lawson {{ rk.label }}
~~~~~~~{% for _ in name %}~{% endfor %}

+ **name:** Lawson {{ rk.label }}
+ **label in ponio:** :cpp:var:`ponio::runge_kutta::l{{ rk.id }}_t`
+ **underlying Runge-Kutta method:** :cpp:type:`ponio::runge_kutta::{{ rk.id }}_t`
+ **stages:** {{ rk.A|length }}
+ **order:** {{ rk.order }}
+ **stages order:** {{ rk.stage_order }}
{% if 'bib' in rk %}+ **bibliography:** [{{ rk.bib.bib }}]({{ rk.bib.url }}){%- endif %}

.. math::

    \begin{array}{c|{%- for ci in rk.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.butcher.A %}
          {{ rk.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.butcher.b|join(' & ') }}
    \end{array}

The stability function:

.. math::

    {{ rk.stability_function }}

{% endif %}{% endfor %}


Embedded methods
----------------

{% for rk in list_erk %}{% if rk.b2 is defined %}

Lawson {{ rk.label }}
~~~~~~~{% for _ in name %}~{% endfor %}

+ **name:** Lawson {{ rk.label }}
+ **label in ponio:** :cpp:var:`ponio::runge_kutta::l{{ rk.id }}_t`**underlying Runge-Kutta method:** :cpp:type:`ponio::runge_kutta::{{ rk.id }}_t`
+ **stages:** {{ rk.A|length }}
+ **order:** {{ rk.order }}
+ **stages order:** {{ rk.stage_order }}
{% if 'bib' in rk %}+ **bibliography:** [{{ rk.bib.bib }}]({{ rk.bib.url }}){%- endif %}

.. math::

    \begin{array}{c|{%- for ci in rk.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.butcher.A %}
          {{ rk.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.butcher.b|join(' & ') }} \\
        \hline
          & {{ rk.butcher.b2|join(' & ') }}
    \end{array}

The stability function:

.. math::

    {{ rk.stability_function }}

{% endif %}{% endfor %}
