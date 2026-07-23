List of additive Runge-Kutta methods
====================================

{% for rk in list_ark %}{% if rk.explicit.b2 is undefined %}

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
          & {{ rk.explicit.butcher.b|join(' & ') }}
    \end{array}

Butcher tableau of implicit part:

.. math::
    \begin{array}{c|{%- for ci in rk.implicit.butcher.c -%}c{%- endfor -%}}
        {%- for ai in rk.implicit.butcher.A %}
          {{ rk.implicit.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
        {%- endfor %}
        \hline
          & {{ rk.implicit.butcher.b|join(' & ') }}
    \end{array}

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

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
