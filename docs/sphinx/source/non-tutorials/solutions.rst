.. _solutions-label:

Solutions to the exercises
**************************

Lennard Jones fluid
===================

Fix a broken input
------------------

.. container:: justify

    You can download the |input_broken_solution| I wrote.
    The key to make the simulation starts, is to initially 
    reduce the *timestep*, as well as the imposed *temperature*.
    Note that in order to make sure that the temperature of the particles
    quickly reaches a low value before the simulation crashes, I also reduced 
    the *damping* parameter of the *Langevin* command:

..  code-block:: lammps

    fix mylgv all langevin 0.001 0.001 0.001 1530917
    timestep 0.0001

.. container:: justify

    Then, after the first *run* finishes, the energy of the system 
    has reduced due to the motion of the atoms, and a second *run*
    with the original *timestep* and *Langevin* parameters can start
    without issue. 

.. container:: justify

    In some case, several consecutive steps can be necessary:

..  code-block:: lammps

    fix mylgv all langevin 0.0001 0.0001 0.001 1530917
    timestep 0.0001
    run 10000

    timestep 0.001
    run 10000

    timestep 0.01
    run 10000

.. container:: justify

    Use trial and error to determine the best approach for
    a given system.

.. |input_broken_solution| raw:: html

    <a href="../../../../inputs/level0/lennard-jones-fluid/exercises/broken/input.lammps" target="_blank">input</a>

Create a demixed dense phase
----------------------------

.. container:: justify

    You can download the |input_demixed_solution| I wrote. Note that 
    I use a large number of particles: 8000 for each type. Then,
    the key to create a demixing phase is to play with the Lennard-Jones 
    parameters:

..  code-block:: lammps

    pair_coeff 1 1 5.0 1.0
    pair_coeff 2 2 5.0 1.0
    pair_coeff 1 2 0.05 1.0

.. container:: justify

    First, notice that both particle types have the same :math:`\sigma` value of 1.0
    so that both particles have the same diameter. Second, note the large energy parameter :math:`\epsilon = 5.0`
    for self interaction (i.e.) interaction between particles of same type, and the low 
    energy parameter :math:`\epsilon = 0.05` for interaction between particles of different types.

.. container:: justify

    Finally, for adjusting the box volume and create a liquid looking phase, the 
    pressure was imposed by replacing *fix nve* by *fix nph*:

..  code-block:: lammps

    fix mynph all nph iso 1.0 1.0 1.0

.. container:: justify

    With *fix nph* and a pressure of 1, LAMMPS adjusts the box dimension until the 
    pressure is close to 1, which is require reducing the initial box dimensions.

.. |input_demixed_solution| raw:: html

    <a href="../../../../inputs/level0/lennard-jones-fluid/exercises/demixion/input.lammps" target="_blank">input</a>

Create dumbbell molecules
-------------------------

.. container:: justify

    You can download the |input_dumbbell_solution| I wrote. The first important 
    change is to choose an *atom_style* that allows for bond creation, and 
    to specify the *bond_style*:

.. |input_dumbbell_solution| raw:: html

    <a href="../../../../inputs/level0/lennard-jones-fluid/exercises/dumbbell/input.lammps" target="_blank">input</a>

..  code-block:: lammps

    atom_style molecular
    bond_style harmonic

.. container:: justify

    When creating the box, it is important to make space in the memory for 
    the bonds:

..  code-block:: lammps

    create_box 2 simulation_box bond/types 2 extra/bond/per/atom 1

.. container:: justify

    Then, import the *molecule templates*, and use these templates
    when creating the atoms:

..  code-block:: lammps

    molecule dumbell1 dumbell1.mol
    molecule dumbell2 dumbell2.mol
    create_atoms 0 random 750 341341 simulation_box mol dumbell1 8766
    create_atoms 0 random 50 678865 simulation_box mol dumbell2 8751

.. container:: justify

    You can download the molecule templates for |mol1_dumbbell_solution|
    and |mol2_dumbbell_solution|. Finally, some parameters for the two
    types of bonds, namely their rigidity and equilibrium lengths is specified:

..  code-block:: lammps

    bond_coeff 1 5 0.5
    bond_coeff 2 5 1.5

.. |mol1_dumbbell_solution| raw:: html

    <a href="../../../../inputs/level0/lennard-jones-fluid/exercises/dumbbell/dumbell1.mol" target="_blank">type-1</a>

.. |mol2_dumbbell_solution| raw:: html

    <a href="../../../../inputs/level0/lennard-jones-fluid/exercises/dumbbell/dumbell2.mol" target="_blank">type-2</a>

Pulling on a carbon nanotube
============================

Plot the strain-stress curves
-----------------------------

.. container:: justify

    You can download the |input_stress_strain_solution1|
    and |input_stress_strain_solution2| I wrote.

.. container:: justify

    The main difficulty here is to calculate the stress. Here, 
    the stress is calculated as the force divided by the 
    surface area of the CNT. Note that the surface area 
    of a CNT is not a well defined quantity. I choose to 
    define the area as the perimeter of the CNT multiplied by the 
    effective width of the carbon atoms. 

.. container:: justify

    Be careful with units, as the force is either in kCal/mol
    or eV, depending of the unit system, respectively *real* or *metal*.

.. |input_stress_strain_solution1| raw:: html

    <a href="../../../../inputs/level1/breaking-a-carbon-nanotube/exercises/stress-strain/breakable-bonds/input.lammps" target="_blank">input</a>

.. |input_stress_strain_solution2| raw:: html

    <a href="../../../../inputs/level1/breaking-a-carbon-nanotube/exercises/stress-strain/unbreakable-bonds/input.lammps" target="_blank">input</a>

Deform a CNT membrane
---------------------

.. container:: justify

    You can download the |input_membrane_solution1| I wrote.
    The CNT can be replicated using the *replicate* command.
    Then, it is important to change the box to triclinic:

..  code-block:: lammps

    change_box all triclinic

.. container:: justify

    Before deforming the system using:

..  code-block:: lammps

    fix muyef all deform 1 xy erate 5e-5

.. |input_membrane_solution1| raw:: html

    <a href="../../../../inputs/level1/breaking-a-carbon-nanotube/exercises/membrane/input.lammps" target="_blank">input</a>

Polymer in water
================

Add salt to the mixture
-----------------------

.. container:: justify

    You can download the |input_PEG_salt|,
    |data_PEG_salt|,
    and |parm_PEG_salt| files I wrote. It is important to 
    make space for the salt by modifying the data file as follow:

..  code-block:: lammps

    (...)
    9 atom types
    (...)

.. container:: justify

    Additional *mass* and *pair_coeff* lines 
    must be added to the parm file (be careful to use the 
    appropriate units):

..  code-block:: lammps

    (...)
    mass 8 22.98 # Na
    mass 9 35.453 # Cl
    (...)
    pair_coeff 8 8 0.04690 2.43 # Na
    pair_coeff 9 9 0.1500 4.045
    (...)

.. container:: justify

    Finally, here I choose to add the ions using two separate
    *create_atoms* commands with a very small *overlap*
    values, followed by an energy minimization. 

.. |input_PEG_salt| raw:: html

    <a href="../../../../inputs/level2/polymer-in-water/exercises/salt/input.lammps" target="_blank">input</a>

.. |data_PEG_salt| raw:: html

    <a href="../../../../inputs/level2/polymer-in-water/exercises/salt/mix-with-salt.data" target="_blank">data</a>

.. |parm_PEG_salt| raw:: html

    <a href="../../../../inputs/level2/polymer-in-water/exercises/salt/PARM-with-salt.lammps" target="_blank">parm</a>

Coming soon
-----------

.. container:: justify

    This page is currently being written, all solutions will appear here eventually.
