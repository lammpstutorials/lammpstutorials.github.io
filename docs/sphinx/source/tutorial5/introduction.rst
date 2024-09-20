.. figure:: avatars/SiO_gif_light.webp
    :height: 250
    :alt: Figure showing silicon dioxide structure with colored charges as simulated with lammps and reaxff
    :class: only-light
    :align: right

.. figure:: avatars/SiO_gif_dark.webp
    :height: 250
    :alt: Figure showing silicon dioxide structure with colored charges as simulated with lammps and reaxff
    :class: only-dark
    :align: right

The objective of this tutorial is to use a 
reactive force field (ReaxFF :cite:`van2001reaxff, zou2012investigation`)
to calculate the partial charges of a system undergoing
deformation, as well as chemical bond formation and breaking.  

The system simulated here is a block of silicon dioxide (:math:`\text{SiO}_2`) that is deformed 
until rupture. Particular attention is given to the evolution of the atomic charges
during the deformation of the structure, and 
the chemical reactions occurring due to the deformation
are tracked over time.
