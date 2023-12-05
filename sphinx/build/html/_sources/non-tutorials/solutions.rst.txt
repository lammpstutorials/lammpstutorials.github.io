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

Extract radial distribution function
------------------------------------

.. container:: justify

    You can download the |input_PEG_RDF| file I wrote. I use 
    the *compute rdf* command of LAMMPS. I define two different *compute*,
    one for the H2O-H2O RDF, and another one for the H2O-PEG RDF:

..  code-block:: lammps
        
    # H2O (type 1) - H2O (type 1)
    compute myRDF_H2O_H2O all rdf 200 1 1 cutoff 10
    fix myat1 all ave/time 10 20000 200000 c_myRDF_H2O_H2O[*] file H2O-H2O.dat mode vector

    # PEG (type 3, 4, and 6) - H2O (type 1)
    compute myRDF_PEG_H2O all rdf 200 3 1 4 1 6 1 cutoff 10
    fix myat2 all ave/time 10 20000 200000 c_myRDF_PEG_H2O[*] file PEG-H2O.dat mode vector
    
.. |input_PEG_RDF| raw:: html

    <a href="../../../../inputs/level2/polymer-in-water/exercises/radial-distribution-function/input.lammps" target="_blank">input</a>

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

Nanosheared electrolyte
=======================

Induce a Poiseuille flow
------------------------

.. container:: justify

    Here the *input* scripts written during the last part *Imposed shearing* of the
    tutorial is adapted so that, instead of a shearing induced by the wall, the fluid is moving
    thanks through an additional force applied to all the water molecules and ions.
    
.. container:: justify

    To do so, here are the most important commands used to properly thermalize the system:

..  code-block:: lammps
        
    fix mynve all nve
    compute tliq fluid temp/partial 0 1 1
    fix myber1 fluid temp/berendsen 300 300 100
    fix_modify myber1 temp tliq
    compute twall wall temp
    fix myber2 wall temp/berendsen 300 300 100
    fix_modify myber2 temp twall
    fix myshk H2O shake 1.0e-4 200 0 b 1 a 1

.. container:: justify

    Note that here, walls wont move and they can be thermalized in all 3 directions.
    There is no need for recentering as the walls here, but instead one can keep them 
    in place by adding springs to every atom of the walls:

..  code-block:: lammps

    fix myspring wall spring/self 10.0 xyz

.. container:: justify

    Finally, let us apply a force to the fluid group along the :math:`x`
    direction:

..  code-block:: lammps

    fix myadf fluid addforce 3e-2 0.0 0.0

.. container:: justify

    One can have a look at the velocity profiles. The fluid shows the characteristic
    parabolic shape of Poiseuille flow in the case of a non-slip solid surface.
    To obtain a a smooth curve, I ran the simulation for a total duration of :math:`1\,\text{ns}` which takes time. 
    You can use a lower duration like :math:`100\,\text{ps}` and still obtain reasonable results.

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/shearing-poiseuille-light.png
    :alt: Velocity of the fluid forming a Poiseuille flow
    :class: only-light

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/shearing-poiseuille-dark.png
    :alt: Velocity of the fluid forming a Poiseuille flow
    :class: only-dark

..  container:: figurelegend

    Figure: Velocity profiles for water molecules and walls along the *z* axis.
    The line is the Stokes equation, and the disk are the results of the simulations.
    
.. container:: justify

    The fitting of the velocity profile was made using:

.. math::

    v = - \alpha \dfrac{F_\text{rho}}{\eta} \left( \dfrac{z^2}{2} - \dfrac{h^2}{8} \right)

.. container:: justify

    Where :math:`F_\text{rho} = f \time \rho` is the force per unit volume of fluid,
    :math:`\eta` the fluid viscosity,
    :math:`h = 1.6\,\text{nm}` the pore size.
    A small correction :math:`\alpha = 0.78` was used. This correction 
    compensates the fact that using bulk density and bulk viscosity is 
    wrong in such nanoconfined pore.
    Note also that the position of the interface, which sets the pore size
    :math:`h = 1.6\,\text{nm}`, is an important parameter that in principle 
    must be chosen with a lot of care, but here goes beyond the scope of this 
    exercise.

.. container:: justify
    
    Here, the force value :math:`f = 0.03\,\text{kcal/mol/Å}` is a reasonable choice
    that was chosen after some prior calibration (see below). :math:`f` is
    small enough so that the system remains in the linear response regime,
    but also large enough so that one can easily differentiate
    the signal from the noise.

.. container:: justify

    The first and most important technical difficulty of any
    out-of-equilibrium simulation is to choose the value of the forcing :math:`f`.
    If the forcing is too large, the system may not be in a linear response regime,
    meaning that the results are forcing-dependent, thus likely quite meaningless. If
    the forcing is too small, the motion of the system will be difficult to measure
    due to the low signal-to-noise ratio.

.. container:: justify

    In the present case, one can perform a calibration by running several simulations 
    with different force values :math:`f`, and then by plotting the velocity of
    the center of mass :math:`v_\text{cm}` of the fluid as a function of the force.
    Here, I present the results I have obtained by performing the simulations with 
    different values of the forcing. :math:`v_\text{cm}` can be extracted by adding the following command
    to the *input*:

..  code-block:: lammps

    variable vcm_fluid equal vcm(fluid,x)
    fix myat1 all ave/time 10 100 1000 v_vcm_fluid file vcm_fluid.dat

.. container:: justify

    The results I have obtained show that as long as the force is lower
    than about :math:`0.04\,\text{kcal/mol/Å}`, there is reasonable linearity
    between force and fluid velocity.

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/calibration-force-light.png
    :alt: Velocity of the fluid under imposed force (POISEUILLE FLOW)
    :class: only-light

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/calibration-force-dark.png
    :alt: Velocity of the fluid under imposed force (POISEUILLE FLOW)
    :class: only-dark

..  container:: figurelegend

    Figure: Ratio between the velocity of the center of mass :math:`v_\text{cm}` of the fluid
    and the forcing :math:`f` as a function of the forcing.
