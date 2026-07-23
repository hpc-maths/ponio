List of embedded Runge-Kutta methods
====================================

Explicit methods
----------------

{% for rk in list_erk %}{% if rk.b2 is defined %}
{{ rk.label }}
~~~~~~~~~~~~~~

+ **name:** {{ rk.label }}
+ **label in ponio:** :cpp:type:`ponio::runge_kutta::{{ rk.id }}_t`
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

Diagonal implicit methods
-------------------------

{% for rk in list_dirk %}{% if rk.b2 is defined %}

{{ rk.label }}
~~~~~~~~~~~~~~

+ **name:** {{ rk.label }}
+ **label in ponio:** :cpp:var:`ponio::runge_kutta::{{ rk.id }}_t`
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

Lawson methods
--------------

{% for rk in list_erk %}{% if rk.b2 is defined %}

Lawson {{ rk.label }}
~~~~~~~~~~~~~~~~~~~~~

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

Additive Runge-Kutta methods
----------------------------

{% for rk in list_ark %}{% if rk.explicit.b2 is defined %}

{{ rk.label }}
~~~~~~~~~~~~~~

+ **name:** {{ rk.label }}
+ **label in ponio:** :cpp:var:`ponio::runge_kutta::{{ rk.id }}_t`
+ **stages:** {{ rk.A|length }}
+ **order:** {{ rk.order }}
{% if 'bib' in rk %}+ **bibliography:** [{{ rk.bib.bib }}]({{ rk.bib.url }}){%- endif %}


Butcher tableau of explicit part:

.. math::
    \begin{array}{c|{%- for ci in rk.explicit.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.explicit.butcher.A %}
          {{ rk.explicit.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.explicit.butcher.b|join(' & ') }} \\
        \hline
          & {{ rk.explicit.butcher.b2|join(' & ') }}
    \end{array}

Butcher tableau of implicit part:

.. math::
    \begin{array}{c|{%- for ci in rk.implicit.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.implicit.butcher.A %}
          {{ rk.implicit.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.implicit.butcher.b|join(' & ') }} \\
        \hline
          & {{ rk.implicit.butcher.b2|join(' & ') }}
    \end{array}

{% endif %}{% endfor %}
