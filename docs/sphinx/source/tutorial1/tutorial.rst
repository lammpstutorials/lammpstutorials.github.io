My first input
==============

To run a simulation using LAMMPS, you need to write an input script
containing a series of commands for LAMMPS to execute, similar to Python
or Bash scripts.  For clarity, the input scripts for this tutorial will
be divided into five categories, which will be filled out step by step.
To set up this tutorial, select *Start LAMMPS Tutorial 1* from
the *Tutorials* menu of LAMMPS--GUI, and follow the
instructions.  This will select (or create, if needed) a folder, place
the initial input file *initial.lmp* in it, and open the file in
the LAMMPS--GUI Editor window.  It should display the following
content:

.. code-block:: lammps

    # PART A - ENERGY MINIMIZATION
    # 1) Initialization
    # 2) System definition
    # 3) Settings
    # 4) Visualization
    # 5) Run

Everything that appears after a hash symbol (#) is a comment
and ignored by LAMMPS.  LAMMPS-GUI will color such comments in red.
These five categories are not required in every input script and do not
necessarily need to be in that exact order.  For instance, the *Settings*
and the *Visualization* categories could be inverted, or
the *Visualization* category could be omitted.  However, note that
LAMMPS reads input files from top to bottom and processes each command
*immediately*.  Therefore, the *Initialization* and
*System definition* categories must appear at the top of the
input, and the *Run* category must appear at the bottom.  Also, the
specifics of some commands can change after global settings are modified, so the
order of commands in the input script is important.

Initialization
--------------

In the first section of the script, called *Initialization*,
global parameters for the simulation are defined, such as units, boundary conditions
(e.g., periodic or non-periodic), and atom types (e.g., uncharged point particles
or extended spheres with a radius and angular velocities). These commands must be
executed *before* creating the simulation box or they will cause
an error. Similarly, many LAMMPS commands may only be
entered *after* the simulation box is defined. Only a limited
number of commands may be used in both cases. Update the initial.lmp file
so that the *Initialization* section appears as follows:

.. code-block:: lammps

    # 1) Initialization
    units lj
    dimension 3
    atom_style atomic
    boundary p p p

The first line, *units lj*, specifies the use of *reduced*  
units, where all quantities are dimensionless.  This unit system is a  
popular choice for simulations that explore general statistical  
mechanical principles, as it emphasizes relative differences between  
parameters rather than representing any specific material.  The second  
line, *dimension 3*, specifies that the simulation is conducted  
in 3D space, as opposed to 2D, where atoms are confined to move only in  
the xy-plane.  The third line, *atom_style atomic*, designates  
the atomic style for representing simple, individual point particles.  
In this style, each particle is treated as a point with a mass, making  
it the most basic atom style.  Other atom styles can incorporate  
additional attributes for atoms, such as charges, bonds, or molecule  
IDs, depending on the requirements of the simulated model.  The last  
line, *boundary p p p*, indicates that periodic boundary  
conditions are applied along all three directions of space, where the  
three p stand for :math:`x`, :math:`y`, and :math:`z`, respectively.  
Alternatives are fixed non-periodic (f), shrink-wrapped non-periodic (s), and  
shrink-wrapped non-periodic with minimum (m).  For non-periodic  
boundaries, different options can be assigned to each dimension, making  
configurations like *boundary p p fm* valid for systems such as  
slab geometries.

.. admonition:: Note
    :class: info

    Strictly speaking, none of the four commands specified in the
    Initialization section are mandatory, as they correspond to the
    default settings for their respective global properties.  However,
    explicitly specifying these defaults is considered good practice to
    avoid confusion when sharing input files with other LAMMPS users.

Each LAMMPS command is accompanied by extensive online documentation  
that details the different options for that command.  From the  
LAMMPS--GUI editor buffer, you can access the documentation by  
right-clicking on a line containing a command (e.g., *units lj*)  
and selecting *View Documentation for `units'*.  This action  
should prompt your web browser to open the corresponding URL for the  
online manual.  A screenshot of this context menu is shown in  
Fig. INCLUDE FIGURE

**************************************

The next step is to create the simulation box and populate it with  
atoms.  Modify the *System definition* category of  
*initial.lmp* as shown below:

.. code-block:: lammps

    # 2) System definition
    region simbox block -20 20 -20 20 -20 20
    create_box 2 simbox
    create_atoms 1 random 1500 34134 simbox overlap 0.3
    create_atoms 2 random 100 12756 simbox overlap 0.3

The first line, *region simbox (...)*, defines a region named  
*simbox* that is a block (i.e.,~a rectangular cuboid) extending  
from -20 to 20 units along all three spatial dimensions.  The second  
line, *create_box 2 simbox*, initializes a simulation box based  
on the region *simbox* and reserves space for two types of atoms.

.. admonition:: Note
    :class: info

    From this point on, any command referencing an atom type larger than 2
    will trigger an error.  While it is possible to allocate more atom
    types than needed, you must assign a mass and provide force field
    parameters for each atom type.  Failing to do so will cause LAMMPS to
    terminate with an error.

The third line, *create_atoms (...)*, generates 1500 atoms of  
type 1 at random positions within the *simbox* region.  The  
integer 34134 is a seed for the internal random number generator, which  
can be changed to produce different sequences of random numbers and,  
consequently, different initial atom positions.  The fourth line adds  
100 atoms of type 2.  Both *create_atoms* commands use the  
optional argument *overlap 0.3*, which enforces a minimum  
distance of 0.3 units between the randomly placed atoms.  This  
constraint helps avoid close contacts between atoms, which can lead  
to excessively large forces and simulation instability.

Settings
--------

Next, we specify the settings for the two atom types.  Modify the
*Settings* category of *initial.lmp* as follows:

.. code-block:: lammps

    # 3) Settings
    mass 1 1.0
    mass 2 5.0
    pair_style lj/cut 4.0
    pair_coeff 1 1 1.0 1.0
    pair_coeff 2 2 0.5 3.0

The two *mass* commands assign a mass of 1.0 and 5.0 units to the
atoms of type 1 and 2, respectively.  The third line,
*pair_style lj/cut 4.0*, specifies that the atoms will be
interacting through a Lennard-Jones (LJ) potential with a cut-off equal
to :math:`r_c = 4.0` length units :cite:`wang2020lennard,fischer2023history`:

.. math::
   E_{ij}(r) = 4 \epsilon_{ij} \left[ \left( \dfrac{\sigma_{ij}}{r} \right)^{12}
   - \left( \dfrac{\sigma_{ij}}{r} \right)^{6} \right], \quad \text{for} \quad r < r_c,

where :math:`r` is the inter-particle distance, :math:`\epsilon_{ij}` is
the depth of the potential well that determines the interaction strength, and
:math:`\sigma_{ij}` is the distance at which the potential energy equals zero.
The indexes :math:`i` and :math:`j` refer to pairs of particle types.
The fourth line, *pair_coeff 1 1 1.0 1.0*, specifies the
Lennard-Jones coefficients for interactions between pairs of atoms
of type 1: the energy parameter :math:`\epsilon_{11} = 1.0` and
the distance parameter :math:`\sigma_{11} = 1.0`.  Similarly, the last line
sets the Lennard-Jones coefficients for interactions between atoms
of type 2, :math:`\epsilon_{22} = 0.5`, and :math:`\sigma_{22} = 3.0`.

.. admonition:: Note
    :class: info

    By default, LAMMPS calculates the cross coefficients for different atom
    types using geometric average: :math:`\epsilon_{ij} = \sqrt{\epsilon_{ii} \epsilon_{jj}}`,
    :math:`\sigma_{ij} = \sqrt{\sigma_{ii} \sigma_{jj}}`.  In the present case,
    :math:`\epsilon_{12} = \sqrt{1.0 \times 0.5} = 0.707`, and
    :math:`\sigma_{12} = \sqrt{1.0 \times 3.0} = 1.732`.

Single-point energy
^^^^^^^^^^^^^^^^^^^

The system is now fully parameterized, and the input is ready to compute
forces.  Let us complete the two remaining categories,
*Visualization* and *Run*, by adding the following lines
to *initial.lmp*:

.. code-block:: lammps

    # 4) Visualization
    thermo 10
    thermo_style custom step etotal press
    # 5) Run
    run 0 post no