.. figure:: avatars/PEG-dark.webp
    :alt: Movie of a peg polymer molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-dark

.. figure:: avatars/PEG-light.webp
    :alt: Movie of a peg polymer molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-light

The goal of this tutorial is to use LAMMPS to solvate a small
hydrophilic polymer (PEG - polyethylene glycol) in a reservoir of water. 

Once the water reservoir is properly equilibrated
at the desired temperature and pressure, the polymer molecule is added
and a constant stretching force is applied to both ends of the polymer.
The evolution of the polymer length is measured as a function of time.
The GROMOS 54A7 force field :cite:`schmid2011definition` is used for the
PEG, the SPC/Fw model :cite:`wu2006flexible` is used for the water, and
the long-range Coulomb interactions are solved using the PPPM
solver :cite:`luty1996calculating`.

This tutorial was inspired by a
publication by Liese and coworkers, in which molecular dynamics
simulations are compared with force spectroscopy experiments, see
Ref. :cite:`liese2017hydration`.

.. admonition:: Note
    :class: non-title-info

    When mixing different force fields, as is done here with GROMOS 
    and SPC/Fw, users should exercise caution.  The choices made in these tutorials 
    prioritize progressive learning of LAMMPS functionality 
    over strict physical accuracy.  While GROMOS is commonly used with water 
    models from the SPC family :cite:`oostenbrink2004biomolecular`, 
    the inter-compatibility of force fields is not generally guaranteed.
