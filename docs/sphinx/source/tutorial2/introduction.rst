.. figure:: avatars/CNT_dark.webp
    :alt: carbon nanotube image in vacuum
    :height: 250
    :align: right
    :class: only-dark

.. figure:: avatars/CNT_light.webp
    :alt: carbon nanotube image in vacuum
    :height: 250
    :align: right
    :class: only-light

In this tutorial, the system of interest is a small, single-walled carbon
nanotube (CNT) in an empty box.  The CNT is strained
by imposing a constant velocity on the edge atoms.  To illustrate the
difference between conventional and reactive force fields, this tutorial
is divided into two parts: in the first part, a conventional molecular
force field (called OPLS-AA :cite:`jorgensenDevelopmentTestingOPLS1996`)
is used and the bonds between the atoms of the CNT are unbreakable.  In
the second part, a reactive force field (called AIREBO :cite:`stuart2000reactive`)
is used, which allows chemical bonds to break under large strain.

To set up this tutorial, select *Start Tutorial 2* from the *Tutorials*
menu of LAMMPS-GUI and follow the instructions. This will select a folder,
create one if necessary, and place several files into it.  The initial
input file, set up for a single-point energy calculation, will also be
loaded into the editor under the name *unbreakable.lmp*.  Additional files
are a data file containing the CNT topology and geometry, named
*unbreakable.data*, a parameters file named *unbreakable.inc*, as well as
the scripts required for the second part of the tutorial.
