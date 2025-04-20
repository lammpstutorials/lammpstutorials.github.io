.. figure:: avatars/avatar.webp
    :height: 250
    :alt: Carbon nanotube (CNT) embedded in a polymer melt composed of polystyrene with LAMMPS
    :class: only-light
    :align: right

.. figure:: avatars/avatar-dm.webp
    :height: 250
    :alt: Carbon nanotube (CNT) embedded in a polymer melt composed of polystyrene with LAMMPS
    :class: only-dark
    :align: right

The goal of this tutorial is to create a model of a carbon nanotube (CNT)
embedded in a polymer melt composed of polystyrene (PS). The REACTER
protocol is used to simulate the polymerization of styrene monomers, and the
polymerization reaction is tracked over time :cite:`gissinger2017polymer,
gissinger2020reacter, gissinger2024molecular`. In contrast to AIREBO
:ref:`carbon-nanotube-label` and ReaxFF :ref:`reactive-silicon-dioxide-label`,
the REACTER protocol relies on the use of a *classical* force field.