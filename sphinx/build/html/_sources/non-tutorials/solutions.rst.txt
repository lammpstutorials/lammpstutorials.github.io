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

    <a href="../../../../inputs/level1/lennard-jones-fluid/exercises/broken/input.lammps" target="_blank">input</a>

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

    <a href="../../../../inputs/level1/lennard-jones-fluid/exercises/demixion/input.lammps" target="_blank">input</a>

Create dumbbell molecules
-------------------------

.. container:: justify

    You can download the |input_dumbbell_solution| I wrote. The first important 
    change is to choose an *atom_style* that allows for bond creation, and 
    to specify the *bond_style*:

.. |input_dumbbell_solution| raw:: html

    <a href="../../../../inputs/level1/lennard-jones-fluid/exercises/dumbbell/input.lammps" target="_blank">input</a>

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

    <a href="../../../../inputs/level1/lennard-jones-fluid/exercises/dumbbell/dumbell1.mol" target="_blank">type-1</a>

.. |mol2_dumbbell_solution| raw:: html

    <a href="../../../../inputs/level1/lennard-jones-fluid/exercises/dumbbell/dumbell2.mol" target="_blank">type-2</a>

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

Make a hydrophobic nanopore
---------------------------

.. container:: justify

    By default, the walls are very hydrophilic due to the 
    *pair_coeff* that are in use:

..  code-block:: lammps

    pair_coeff 1 1 0.185199 3.1589 # water
    (...)
    pair_coeff 5 5 11.697 2.574 # wall

.. container:: justify

    LAMMPS calculate the cross-coefficient between water and wall 
    using geometric average, so
    :math:`\epsilon_\text{1-5} = (0.185199+11.697)/2 = 5.941\,\text{kcal/mol}` and
    :math:`\sigma_\text{1-5} = (3.1589+2.574)/2 = 2.866\,\text{Å}`.

.. container:: justify

    The value :math:`\epsilon_\text{1-5} = 5.941\,\text{kcal/mol}` is extremely high 
    (compare to the water-water energy :math:`\epsilon_\text{1-1} = 0.185199\,\text{kcal/mol}`),
    which is the reason for the extremely hydrophilic surfaces.

.. container:: justify

    The walls can be made more hydrophobic by reducing the 
    LJ energy of interaction :math:`\epsilon_\text{1-5}`.
    This can be done by writing a new line for the *pair_coeff* between
    atoms of type 1 (water oxygen) and wall (type 5):

..  code-block:: lammps

    pair_coeff 1 5 0.0059411 2.86645

.. container:: justify

    Here I choose to reduce :math:`\epsilon_\text{1-5}` by a factor of 1000,
    which makes it lower than the water-water energy :math:`\epsilon_\text{1-1}`,
    thus making the walls more hydrophobic.

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/hydrophobic-pore-light.png
    :alt: hydrophobic vs hydrophilic pores : density profiles
    :class: only-light

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/hydrophobic-pore-dark.png
    :alt: hydrophobic vs hydrophilic pores : density profiles
    :class: only-dark

..  container:: figurelegend

    Figure: Density profile for the water along the *z* axis
    comparing the original hydrophilic pore with the hydrophobic pore 
    with :math:`\epsilon_\text{1-5} = 0.005941\,\text{kcal/mol}`.

Induce a Poiseuille flow
------------------------

.. container:: justify

    Here the *input* script written during the last part *Imposed shearing* of the
    tutorial is adapted so that, instead of a shearing induced by the relative motion of the walls,
    the fluid motion is generated by an additional force applied to both water molecules and ions.
    
.. container:: justify

    To do so, here are the most important commands used to properly
    thermalize the system:

..  code-block:: lammps
        
    fix mynve all nve
    compute tliq fluid temp/partial 0 1 1
    fix myber1 fluid temp/berendsen 300 300 100
    fix_modify myber1 temp tliq
    compute twall wall temp
    fix myber2 wall temp/berendsen 300 300 100
    fix_modify myber2 temp twall

.. container:: justify

    Here, since walls wont move, they can be thermalized in all
    3 directions and there is
    no need for recentering. Instead, one can keep the walls 
    in place by adding springs to every atom:

..  code-block:: lammps

    fix myspring wall spring/self 10.0 xyz

.. container:: justify

    Finally, let us apply a force to the fluid group along the :math:`x`
    direction:

..  code-block:: lammps

    fix myadf fluid addforce 3e-2 0.0 0.0

.. container:: justify

    The choice of a force equal to :math:`f = 0.03\,\text{kcal/mol/Å}`
    is discussed below.

.. container:: justify

    One can have a look at the velocity profiles. The fluid shows the characteristic
    parabolic shape of Poiseuille flow in the case of a non-slip solid surface.
    To obtain smooth looking data, I ran the simulation for a total duration of :math:`1\,\text{ns}`. 
    To lower the duration of the computation, don't hesitate to
    use a shorter duration like :math:`100\,\text{ps}`.

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/shearing-poiseuille-light.png
    :alt: Velocity of the fluid forming a Poiseuille flow
    :class: only-light

.. figure:: ../tutorials/figures/level2/nanosheared-electrolyte/shearing-poiseuille-dark.png
    :alt: Velocity of the fluid forming a Poiseuille flow
    :class: only-dark

..  container:: figurelegend

    Figure: Velocity profiles of the water molecules along the *z* axis (disks).
    The line is the Poiseuille equation.
    
.. container:: justify

    The fitting of the velocity profile was made using the following Poiseuille equation,

.. math::

    v = - \alpha \dfrac{f \rho}{\eta} \left( \dfrac{z^2}{2} - \dfrac{h^2}{8} \right),

.. container:: justify

    which can be derived from the Stokes equation :math:`\eta \nabla \textbf{v} = - \textbf{f} \rho`
    where :math:`f` is the applied force,
    :math:`\rho` is the fluid density,
    :math:`\eta` is the fluid viscosity, and
    :math:`h = 1.2\,\text{nm}` is the pore size.
    A small correction :math:`\alpha = 0.78` was used. This correction 
    compensates the fact that using bulk density and bulk viscosity is obviously
    no correct in such nanoconfined pore. More subtle corrections could be applied
    by correcting both density and viscosity based on independent measurement, but this is 
    beyond the scope of the present exercise.

.. container:: justify

    **Choosing the right force**

.. container:: justify

    The first and most important technical difficulty of any
    out-of-equilibrium simulation is to choose the value of the forcing :math:`f`.
    If the forcing is too large, the system may not be in a linear response regime,
    meaning that the results are forcing-dependent (and likely quite meaningless). If
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

    The results  show that as long as the force is lower
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
    and the forcing :math:`f` as a function of the forcing


Water adsorption in silica
==========================

Mixture adsorption
------------------

.. container:: justify

    You can download the |input_mixture| for the combine water and CO2
    adsorption.
    One of the first step is to create both type of molecules
    before starting the GCMC:

..  code-block:: lammps

    molecule h2omol H2O.mol
    molecule co2mol CO2.mol
    create_atoms 0 random 5 456415 NULL mol h2omol 454756 overlap 2.0 maxtry 50
    create_atoms 0 random 5 373823 NULL mol co2mol 989812 overlap 2.0 maxtry 50

.. container:: justify

    One must be careful to properly write the parameters of the system,
    with all the proper cross coefficients:

..  code-block:: lammps

    pair_coeff * * vashishta ../../Potential/SiO.1990.vashishta Si O NULL NULL NULL NULL
    pair_coeff * * lj/cut/tip4p/long 0 0
    pair_coeff 1 3 lj/cut/tip4p/long 0.0057 4.42 # epsilonSi = 0.00403, sigmaSi = 3.69
    pair_coeff 1 5 lj/cut/tip4p/long 0.01096 3.158 # epsilonSi = 0.00403, sigmaSi = 3.69
    pair_coeff 1 6 lj/cut/tip4p/long 0.007315 3.2507 # epsilonSi = 0.00403, sigmaSi = 3.69
    pair_coeff 2 3 lj/cut/tip4p/long 0.0043 3.12 # epsilonO = 0.0023, sigmaO = 3.091
    pair_coeff 2 5 lj/cut/tip4p/long 0.0101 2.858 # epsilonO = 0.0023, sigmaO = 3.091
    pair_coeff 2 6 lj/cut/tip4p/long 0.0065 2.9512 # epsilonO = 0.0023, sigmaO = 3.091
    pair_coeff 3 3 lj/cut/tip4p/long 0.008 3.1589
    pair_coeff 3 5 lj/cut/tip4p/long 0.01295 2.8924
    pair_coeff 3 6 lj/cut/tip4p/long 0.0093 2.985
    pair_coeff 4 4 lj/cut/tip4p/long 0.0 0.0
    pair_coeff 5 5 lj/cut/tip4p/long 0.0179 2.625854
    pair_coeff 6 6 lj/cut/tip4p/long 0.0106 2.8114421 

.. container:: justify

    Here, I choose to thermalize all species separately:

..  code-block:: lammps

    compute ctH2O H2O temp
    compute_modify ctH2O dynamic yes
    fix mynvt1 H2O nvt temp 300 300 0.1
    fix_modify mynvt1 temp ctH2O

    compute ctCO2 CO2 temp
    compute_modify ctCO2 dynamic yes
    fix mynvt2 CO2 nvt temp 300 300 0.1
    fix_modify mynvt2 temp ctCO2

    compute ctSiO SiO temp
    fix mynvt3 SiO nvt temp 300 300 0.1
    fix_modify mynvt3 temp ctSiO

.. container:: justify

    Finally, adsorption is made with two separates *fix gcmc* commands
    placed in a loop: 

..  code-block:: lammps

    label loop
    variable a loop 30

    fix fgcmc_H2O H2O gcmc 100 100 0 0 65899 300 -0.5 0.1 mol h2omol tfac_insert ${tfac} group H2O shake shak full_energy pressure 100 region system
    run 500
    unfix fgcmc_H2O

    fix fgcmc_CO2 CO2 gcmc 100 100 0 0 87787 300 -0.5 0.1 mol co2mol tfac_insert ${tfac} group CO2 full_energy pressure 100 region system
    run 500
    unfix fgcmc_CO2

    next a
    jump SELF loop

.. container:: justify

    Here I choose to apply the first *fix gcmc* for the :math:`\text{H}_2\text{O}` for 500 steps,
    then unfix it before starting the second *fix gcmc* for the :math:`\text{CO}_2` for 500 steps as well.
    Then, thanks to the *jump*, these two fixes are applied successively 30 times each, allowing for the 
    progressive adsorption of both species.

.. |input_mixture| raw:: html

    <a href="../../../../inputs/level3/water-adsorption-in-silica/Exercises/MixtureH2OCO2/input.lammps" target="_blank">input</a>

Adsorb water in ZIF-8 nanopores
-------------------------------

.. container:: justify

    You can download the |input_zif| for the water adsorption in Zif-8,
    which you have to place in the same folder as the *zif-8.data*,
    *parm.lammps*,
    and *water.mol* files.

.. |input_zif| raw:: html

    <a href="../../../../inputs/level3/water-adsorption-in-silica/Exercises/Zif-8/input.lammps" target="_blank">input</a>

.. container:: justify

    Apart from the parameters and topology, the *input* is
    quite similar to the one developped in the case of the crack
    silica.

.. container:: justify

    You should observe an increase of the number of molecule with time.
    Run much longer simulation if you want to saturate the porous material
    with water.

.. figure:: ../tutorials/figures/level3/water-adsorption-in-silica/number_evolution_zif-light.png
    :alt: Water molecule in Zif material with GCMC in LAMMPS
    :class: only-light

.. figure:: ../tutorials/figures/level3/water-adsorption-in-silica/number_evolution_zif-dark.png
    :alt: Water molecule in Zif material with GCMC in LAMMPS
    :class: only-dark

..  container:: figurelegend

    Figure: Number of water molecule in Zif-8 during the first :math:`10\,ps`.

Reactive silicon dioxide
========================

Add O2 molecules
----------------

.. container:: justify

    In a separate folder, create a new input file,
    and copy the same first lines as previously in it
    (just adapt the path to *silica-deformed.data* accordingly): 

..  code-block:: lammps

    units real
    atom_style full

    read_data ../../Deform/silica-deformed.data

    mass 1 28.0855 # Si
    mass 2 15.999 # O

    pair_style reaxff NULL safezone 3.0 mincap 150
    pair_coeff * * ../RelaxSilica/reaxCHOFe.ff Si O
    fix myqeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff maxiter 400

..  container:: justify

    Optionally, let us shift the structure to recenter it in the box. The best value 
    for the shift may be different in your case. This step is not necessary, but the
    recentered system looks better.

..  code-block:: lammps

    displace_atoms all move -13 0 0 units box

..  container:: justify

    Then, let us import the molecule template *O2.mol* and create 10 molecules. 
    The overlap and maxtry keywords allow us to prevent overlapping
    between the atoms:

..  code-block:: lammps

    molecule O2mol O2.mol
    create_atoms 0 random 10 456415 NULL mol O2mol 454756 overlap 3.0 maxtry 50

..  container:: justify

    The value of 3 Angstroms for the minimum interatomic overlapping is 
    very safe for the present system. Smaller values may lead to molecules being 
    too close from each others.

..  container:: justify

    Finally, let us minimize the energy of the system, and run for :math:`10\,\text{ps}`:

..  code-block:: lammps

    minimize 1.0e-4 1.0e-6 100 1000
    reset_timestep 0

    group grpSi type 1
    group grpO type 2
    variable totqSi equal charge(grpSi)
    variable totqO equal charge(grpO)
    variable nSi equal count(grpSi)
    variable nO equal count(grpO)
    variable qSi equal v_totqSi/${nSi}
    variable qO equal v_totqO/${nO}

    dump dmp all custom 100 dump.lammpstrj id type q x y z
    thermo 5
    thermo_style custom step temp etotal press vol v_qSi v_qO
    fix myspec all reaxff/species 5 1 5 species.log element Si O

    fix mynvt all nvt temp 300.0 300.0 100
    timestep 0.5 

    run 20000

..  container:: justify

    You can vizualise the :math:`\text{O}_2` molecules with VMD, or have a look at the
    *species.log* file:

..  code-block:: lammps

    #  Timestep    No_Moles    No_Specs   Si192O384          O2
              5          11           2           1          10

..  container:: justify

    One can see that some reactions occur in the system, and
    that eventually some of
    the :math:`\text{O}_2` molecules react and reabsorb on the 
    main structure:

..  code-block:: lammps

    #  Timestep    No_Moles    No_Specs   Si192O388          O2
          20000           9           2           1           8

.. figure:: ../tutorials/figures/level3/reactive-silicon-dioxide/O2_light.png
    :alt: Silicon oxide with additional O2 molecules
    :class: only-light

.. figure:: ../tutorials/figures/level3/reactive-silicon-dioxide/O2_dark.png
    :alt: Silicon oxide with additional O2 molecules
    :class: only-dark

..  container:: figurelegend

    Figure: Deformed structure with some :math:`\text{O}_2` molecules

Decorate dandling oxygens
-------------------------

..  container:: justify

    Space must be made for the hydrogen atoms. Modify the *silica-deformed.data* file
    so that it starts with:

..  code-block:: lammps

    576 atoms
    3 atom types

..  container:: justify

    Also add the mass of the hydrogen:

..  code-block:: lammps

    Masses

    1 28.0855
    2 15.999
    3 1.008

..  container:: justify

    It is also important to change the *pair_coeff*:

..  code-block:: lammps

    pair_coeff * * ../../RelaxSilica/reaxCHOFe.ff Si O H

..  container:: justify

    One can create randomly a few hydrogen atoms:

..  code-block:: lammps

    create_atoms 3 random 10 456415 NULL overlap 3.0 maxtry 50

..  container:: justify

    Equilibrate the system, you should see the hydrogen atoms 
    progressively decorating the surface of the SiO2 structure:

..  code-block:: lammps

    #  Timestep    No_Moles    No_Specs    Si192O384        H
              5          11           2            1       10
    (...)
    #  Timestep    No_Moles    No_Specs Si192O384H10
           5000           1           1            1

.. figure:: ../tutorials/figures/level3/reactive-silicon-dioxide/exercice-light.png
    :alt: Silicon oxide decorated with hydrogens
    :class: only-light

.. figure:: ../tutorials/figures/level3/reactive-silicon-dioxide/exercice-light.png
    :alt: Silicon oxide decorated with hydrogens
    :class: only-dark

..  container:: figurelegend

    Figure: Hydrogen atoms are in white, oxygen in red, and silicon in yellow.