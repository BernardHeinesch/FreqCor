{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :noindex:

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
      :nosignatures:

   {% for item in methods %}
      {{ item }}
   {% endfor %}
   {% endif %}

   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree:
      :nosignatures:

   {% for item in attributes %}
      {{ item }}
   {% endfor %}
   {% endif %}
