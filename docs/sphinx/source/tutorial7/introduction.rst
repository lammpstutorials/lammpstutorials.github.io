.. figure:: figures/avatar_light.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-light
    :align: right

.. figure:: figures/avatar_dark.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-dark
    :align: right

The objective of this tutorial is to measure a free
energy profile of particles across a barrier potential
using two methods: free sampling
and umbrella sampling :cite:`kastner2011umbrella, allen2017computer, frenkel2023understanding`.

For simplicity and to reduce computation time, the barrier potential will
be imposed on the atoms with an additional force, mimicking the presence of
a repulsive area in the middle of the simulation box without the need to
simulate extra atoms. The procedure is valid for more complex systems and
can be adapted to many other situations, such as measuring the adsorption
barrier near an interface, or for calculating a translocation
barrier through a membrane :cite:`wilson1997adsorption, makarov2009computer, gravelle2021adsorption`.
