.. figure:: avatars/nanoconfined-electrolyte-dark.png
    :height: 250
    :alt: Electrolyte nano-confined in a slit pore
    :class: only-dark
    :align: right

.. figure:: avatars/nanoconfined-electrolyte-light.png
    :height: 250
    :alt: Electrolyte nano-confined in a slit pore
    :class: only-light
    :align: right

The objective of this tutorial is to simulate an electrolyte
nanoconfined and sheared between two walls.  The density
and velocity profiles of the fluid in the direction normal to the walls are
extracted to highlight the effect of confining a fluid on its local properties.
This tutorial demonstrates key concepts of combining a fluid and a solid in
the same simulation.  A major difference from the previous tutorial,
:ref:`all-atoms-label`, is that here a rigid four-point
water model named TIP4P/2005 is used :cite:`abascal2005general`.

.. admonition:: Note
    :class: non-title-info
        
    Four-point water models such as TIP4P/2005 are widely used as they offer a
    good compromise between accuracy and computational cost :cite:`kadaoluwa2021systematic`.
