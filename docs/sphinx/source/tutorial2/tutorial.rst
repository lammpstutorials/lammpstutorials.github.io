Unbreakable bonds
=================

With most conventional molecular force fields, the chemical bonds between
atoms are defined at the start of the simulation and remain fixed, regardless
of the forces applied to the atoms.  These bonds are typically modeled as springs
with equilibrium distances :math:`r_0` and force constants :math:`k_\text{b}`:
:math:`U_\text{b} = k_\text{b} \left( r - r_0 \right)^2`.  Additionally, angular and
dihedral constraints are often imposed to preserve the molecular structure
by maintaining the relative orientations of neighboring atoms.

The LAMMPS input
----------------

After completing the setup, the editor should display the following content:

.. code-block:: lammps

    units real
    atom_style molecular
    boundary f f f

    pair_style lj/cut 14.0
    bond_style harmonic
    angle_style harmonic
    dihedral_style opls
    improper_style harmonic
    special_bonds lj 0.0 0.0 0.5

    read_data unbreakable.data
    include unbreakable.inc

    run 0 post no

The chosen unit system is *real* (therefore distances are in
Ångströms (Å), times in femtoseconds (fs), and energies in kcal/mol), the
*atom_style* is *molecular* (therefore atoms are point
particles that can form bonds with each other), and the boundary
conditions are fixed.  The boundary conditions do not matter here, as
the box boundaries were placed far from the CNT.  Just like in the
previous tutorial, :ref:`lennard-jones-label`,
the pair style is *lj/cut* (i.e. a Lennard-Jones potential with
cutoff) and its cutoff is set to 14 Å, which means that only the
atoms closer than this distance interact through the Lennard-Jones
potential.

The *bond_style*, *angle_style*,
*dihedral_style*, and *improper_style* commands specify
the different potentials used to constrain the relative positions of the
atoms.  The *special_bonds* command sets the weighting factors
for the Lennard-Jones interactions between atoms directly connected by
one bond, two bonds, and three bonds, respectively.  This is done for
convenience when parameterizing the force constants for bonds, angles, and
so on.  By excluding the non-bonded (Lennard-Jones) interactions for
these pairs, those interactions do not need to be considered when determining
the force constants.

