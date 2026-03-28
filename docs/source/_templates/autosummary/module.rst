{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}
   :noindex:

{% if functions %}
.. rubric:: Functions

.. autosummary::
   :toctree:
   :nosignatures:

{% for item in functions %}
   {{ item }}
{% endfor %}
{% endif %}

{% if classes %}
.. rubric:: Classes

.. autosummary::
   :toctree:
   :nosignatures:

{% for item in classes %}
   {{ item }}
{% endfor %}
{% endif %}

{% if exceptions %}
.. rubric:: Exceptions

.. autosummary::
   :toctree:

{% for item in exceptions %}
   {{ item }}
{% endfor %}
{% endif %}
