.. figure:: avatars/avatar_light.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-light
    :align: right

.. figure:: avatars/avatar_dark.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-dark
    :align: right

The goal of this tutorial is to create a model of a carbon nanotube (CNT) embedded
in a polymer melt made of polystyrene (PS) (Fig.~\ref{fig:REACT}).  The
REACTER protocol is used to simulate the polymerization of styrene monomers, and the
polymerization reaction is followed in time :cite:`gissinger2017polymer, gissinger2020reacter, gissinger2024molecular`.
In contrast with AIREBO :ref:`carbon-nanotube-label`
and ReaxFF :ref:`reactive-silicon-dioxide-label`, the REACTER
protocol relies on the use of a *classical* force field.