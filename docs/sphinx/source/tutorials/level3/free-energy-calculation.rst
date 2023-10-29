.. _umbrella-sampling-label:

Free energy calculation
***********************

.. container:: hatnote

    Simple sampling of a free energy barrier using umbrella sampling

.. figure:: ../figures/level3/free-energy-calculation/avatar_light.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-light
    :align: right

.. figure:: ../figures/level3/free-energy-calculation/avatar_dark.webp
    :height: 250
    :alt: Lennard Jones atoms simulated with LAMMPS
    :class: only-dark
    :align: right

..  container:: justify

    The objective of this tutorial is to measure the free
    energy profile across a barrier potential using two methods;
    :ref:`free sampling <method1>` and :ref:`umbrella sampling <method2>`.
    
    For the sake of simplicity and in order to reduce the computation time, the
    barrier potential will be imposed artificially to the atoms.
    The procedure is valid for more complex
    systems, and can be adapted to many other situations, for instance 
    for measuring adsorption barrier near a wall, or for calculating translocation
    barrier through a membrane.

.. include:: ../../contact/recommand-lj.rst

.. include:: ../../contact/needhelp.rst

.. include:: ../../contact/2Aug2023.rst

.. _method1:

Method 1: Free sampling
=======================

..  container:: justify

    The most direct way to calculate a free energy profile is to extract
    the partition function from a classic (unbiased) molecular
    dynamics simulation, and then to estimate the Gibbs free
    energy using 
    
.. math:: \Delta G = -RT \ln(p/p_0),
    
..  container:: justify

    where :math:`\Delta G` is the free energy difference, R the
    gas constant, T the temperature, p the
    pressure, and :math:`p_0` the reference pressure.
    As an illustration, let us apply this method to an
    extremely simple configuration that consists in a few
    particles diffusing in a box in presence of a
    position-dependent repealing force that makes the centre
    of the box a relatively unfavourable area to explore.

Basic LAMMPS parameters
-----------------------

..  container:: justify

    Create a folder named *FreeSampling/*, and create an input script
    named *input.lammps* in it. Copy the following lines:

..  code-block:: lammps
    :caption: *to be copied in FreeSampling/input.lammps*

    # define some variables
    variable sigma equal 3.405 # Angstrom
    variable epsilon equal 0.238 # Kcal/mol
    variable U0 equal 1.5*${epsilon} # Kcal/mol
    variable dlt equal 1.0 # Angstrom
    variable x0 equal 10.0  # Angstrom

    # initialise the simulation
    units real
    atom_style atomic
    pair_style lj/cut 3.822 # 2^(1/6) * 3.405 WCA potential
    pair_modify shift yes
    boundary p p p

..  container:: justify

    Here, we start by defining variables for the Lennard-Jones
    interaction :math:`\sigma` and :math:`\epsilon` and for
    the repulsive potential :math:`U (x)`: :math:`U_0`, :math:`\delta`, and :math:`x_0`, 
    see the analytical expression below.

    The system of unit 'real' (for which energy is in kcal/mol, distance in Ångstrom,
    time in femtosecond) has been chosen for practical reason,
    as the WHAM algorithm we are going to use in the second
    part of the tutorial automatically assumes the energy to
    be in kcal/mol. Atoms will interact through a
    Lennard-Jones potential with a cut-off equal to 
    :math:`\sigma \times 2 ^ {1/6}` (i.e. a WCA repulsive
    potential). The potential is shifted to be equal to 0 at
    the cut-off using the pair_modify.

System creation and settings
----------------------------

..  container:: justify

    Let us define the simulation block and randomly add atoms:

..  code-block:: lammps
    :caption: *to be copied in FreeSampling/input.lammps*

    # define the system
    region myreg block -25 25 -5 5 -25 25
    create_box 1 myreg
    create_atoms 1 random 60 341341 myreg overlap 1.0 maxtry 50

    # settings
    mass * 39.95
    pair_coeff * * ${epsilon} ${sigma}
    neigh_modify every 1 delay 4 check yes

..  container:: justify

    Here I am using the argon's values of the Lennard-Jones parameters :math:`\sigma` and
    :math:`\epsilon`, as well as the mass :math:`m = 39.95`
    grams/mole. 
    
    In the previous subsection, the variables :math:`U_0`, :math:`\delta`, and
    :math:`x_0` were defined. They are used to create the repulsive potential
    restricting the atoms to explore the center of the box: 

.. math::

    U(x) = U_0 \left[ \arctan \left( \dfrac{x+x_0}{\delta} \right) - \arctan \left(\dfrac{x-x_0}{\delta} \right) \right]. 
    
..  container:: justify

    From the derivative of the
    potential with respect to :math:`x`, we obtain the expression
    for the force that will be imposed to the atoms:

.. math::

    F(x)= \dfrac{U_0}{\delta} \left[ \dfrac{1}{(x-x_0)^2/\delta^2+1} - \dfrac{1}{(x+x_0)^2/\delta^2+1} \right].

..  container:: justify

    The potential and force along the :math:`x`
    axis resemble:

.. figure:: ../figures/level3/free-energy-calculation/potential-light.png
   :alt: Imposed potential
   :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/potential-dark.png
   :alt: Averaged density profile
   :class: only-dark

   Potential :math:`U (x)` (top) and force :math:`F (x)` (bottom) imposed to the atoms.

**Energy minimization and equilibration**

..  container:: justify

    Let us apply energy minimization to the system, and then impose
    the force \\(F(x)\\) to all of the atoms in the simulation using the 'addforce' command:

..  code-block:: lammps
    :caption: *to be copied in FreeSampling/input.lammps*

    # --------------------- Run
    minimize 1e-4 1e-6 100 1000
    reset_timestep 0

    variable U atom ${U0}*atan((x+${x0})/${dlt})-${U0}*atan((x-${x0})/${dlt})
    variable F atom ${U0}/((x-${x0})^2/${dlt}^2+1)/${dlt}-${U0}/((x+${x0})^2/${dlt}^2+1)/${dlt}
    fix myadf all addforce v_F 0.0 0.0 energy v_U

..  container:: justify

    Finally, let us combine the fix nve with a Langevin
    thermostat to run a molecular dynamics simulation. With
    these two commands, the MD simulation is effectively in the
    NVT ensemble: constant number of atoms :math:`N`, constant
    volume :math:`V`, and constant temperature :math:`T`. Let us
    perform an equilibration of 500000 steps in total,
    using a timestep of 2 ps (i.e. a total duration of 1
    nanoseconds). To make sure that 1 ns is long enough, let us
    record the evolution of the number of atoms in the central
    (energetically unfavorable) region called *mymes*:

..  code-block:: lammps
    :caption: *to be copied in FreeSampling/input.lammps*

    fix mynve all nve
    fix mylgv all langevin 119.8 119.8 50 1530917

    region mymes block -${x0} ${x0} INF INF INF INF 
    variable n_center equal count(all,mymes)
    fix myat all ave/time 10 50 500 v_n_center file density_evolution.dat

    timestep 2.0
    thermo 10000
    run 500000

Run and data acquisition
------------------------

..  container:: justify

    Finally, let us record the density profile of the atoms
    along the :math:`x` axis using the 'ave/chunk' command. A
    total of ten density profiles will be printed. Step counts are
    reset to 0 to synchronize with the output times of
    density/number, and the fix 'myat' is canceled (it has to be
    canceled before a reset time).

..  code-block:: lammps
    :caption: *to be copied in FreeSampling/input.lammps*

    unfix myat
    reset_timestep 0

    compute cc1 all chunk/atom bin/1d x 0.0 1.0
    fix myac all ave/chunk 10 400000 4000000 cc1 density/number file density_profile_8ns.dat
    dump mydmp all atom 200000 dump.lammpstrj

    thermo 100000
    run 4000000

..  container:: justify

    This simulation with a duration of 8 ns needs a few
    minutes to complete. Feel free to increase the 
    duration of the run for smoother results.
    
    You can visualize the dump file using VMD:

.. figure:: ../figures/level3/free-energy-calculation/system-light.png
   :alt: Lennard jones atoms simulated with LAMMPS MD code
   :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/system-dark.png
   :alt: Lennard jones atoms simulated with LAMMPS MD code
   :class: only-dark

   Notice that the density of atoms is lower in the central part of the box, 
   due to the additional force :math:`F (x)`.

Data analysis
--------------

..  container:: justify

    First, let us make sure that the equilibration duration of 1
    ns is long enough by looking at the 'density_evolution.dat' file:

.. figure:: ../figures/level3/free-energy-calculation/density_evolution-light.png
   :alt: Number of particle in the central region as a function of time
   :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/density_evolution-dark.png
   :alt: Number of particle in the central region as a function of time
   :class: only-dark

   Evolution of the number of atoms in the central region during equilibration. 
   
..  container:: justify

    Here, we can clearly see that the number of atoms in the
    central region, :math:`n_\mathrm{central}`, evolves to its equilibrium value
    after about 0.1 ns.
    
    Let us also plot the equilibrium density profile :math:`\rho`:

.. figure:: ../figures/level3/free-energy-calculation/density_profile-light.png
   :alt: Averaged density profile
   :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/density_profile-dark.png
   :alt: Averaged density profile
   :class: only-dark

   Averaged density profiles for the 8 ns run. 
   The value for the reference density :math:`\rho_\mathrm{bulk} = 0.0033`
   was estimated from the raw density profiles.

..  container:: justify

    Then, let us plot :math:`-R T \ln(\rho/\rho_\mathrm{bulk})` and compare it
    with the imposed (reference) potential :math:`U`:

.. figure:: ../figures/level3/free-energy-calculation/freesampling-potential-light.png
   :alt: Averaged density profile
   :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/freesampling-potential-dark.png
   :alt: Averaged density profile
   :class: only-dark

   Calculated potential :math:`-R T \ln(\rho/\rho_\mathrm{bulk})` compared to imposed potential.
   The calculated potential is in blue.

..  container:: justify

    The agreement with the expected energy profile is reasonable,
    despite some noise in the central part. 

The limits of free sampling
---------------------------

..  container:: justify

    If we increase the value of :math:`U_0`, the average number of
    atoms in the central region will decrease, making it
    difficult to obtain a good resolution for the free energy
    profile.

    In that case, it is better to use the umbrella sampling method
    to extract free energy profiles, see the next section.

.. include:: ../../contact/supportme.rst

.. _method2:

Method 2: Umbrella sampling
===========================

..  container:: justify

    Umbrella sampling is a biased molecular dynamics method,
    i.e. a method in which additional forces are added to the
    atoms in order to make the unfavourable states more likely
    to occur.
    
    Keeping the present configuration, we are going to force one of the atom to
    explore the central region of the box. To do so, we
    are going to add a potential :math:`V` to one
    of the particle, and force it to move along the axe :math:`x`.
    The chosen path is called the axe of reaction. The final
    simulation will be analyzed using the weighted histogram
    analysis method (WHAM), which allows to remove the effect of
    the bias and eventually deduce the unbiased free energy profile.

LAMMPS input script
-------------------

..  container:: justify

    Create a new folder called *BiasedSampling/*, create a new input file 
    named *input.lammps* in it, and copy the following lines:

..  code-block:: lammps
    :caption: *to be copied in BiasedSampling/input.lammps*

    # define a bunch of variables
    variable sigma equal 3.405 # Angstrom
    variable epsilon equal 0.238 # Kcal/mol
    variable U0 equal 10*${epsilon} # Kcal/mol
    variable dlt equal 0.5 # Angstrom
    variable x0 equal 5.0  # Angstrom
    variable k equal 1.5 # Kcal/mol/Angstrom^2

    # initialise the simulation
    units real
    atom_style atomic
    pair_style lj/cut 3.822 # 2^(1/6) * 3.405 WCA potential
    pair_modify shift yes
    boundary p p p

    # define the system
    region myreg block -25 25 -5 5 -25 25
    create_box 2 myreg
    create_atoms 2 single 0 0 0
    create_atoms 1 random 5 341341 myreg

    # settings
    mass * 39.948
    pair_coeff * * ${epsilon} ${sigma}
    neigh_modify every 1 delay 4 check yes
    group topull type 2

    # run
    variable U atom ${U0}*atan((x+${x0})/${dlt})-${U0}*atan((x-${x0})/${dlt})
    variable F atom ${U0}/((x-${x0})^2/${dlt}^2+1)/${dlt}-${U0}/((x+${x0})^2/${dlt}^2+1)/${dlt}
    fix pot all addforce v_F 0.0 0.0 energy v_U

    fix mynve all nve
    fix mylgv all langevin 119.8 119.8 50 1530917
    timestep 2.0
    thermo 100000
    run 500000
    reset_timestep 0

    dump mydmp all atom 1000000 dump.lammpstrj

..  container:: justify

    So far, this code resembles the one of Method 1,
    except for the additional particle of type 2. This
    particle is identical to the particles of type 1 (same
    mass and Lennard-Jones parameters), but will be exposed to the
    biasing potential.

    The value of the potential :math:`U_0` was chosen to be much larger than in part 1, 
    just to proof that umbrella sampling can easily deal with huge potential value,
    while free sampling couldn't.

    Let us create a loop with 50 steps, and move progressively
    the centre of the bias potential by increment of 0.1 nm:

..  code-block:: lammps
    :caption: *to be copied in BiasedSampling/input.lammps*

    variable a loop 50
    label loop
    variable xdes equal ${a}-25
    variable xave equal xcm(topull,x)
    fix mytth topull spring tether ${k} ${xdes} 0 0 0
    run 200000
    fix myat1 all ave/time 10 10 100 v_xave v_xdes file data-k1.5/position.${a}.dat
    run 1000000
    unfix myat1
    next a
    jump SELF loop

..  container:: justify

    A folder named *data-k1.5/* needs to be created within *BiasedSampling/*.

    The spring command serves to impose the
    additional harmonic potential with spring constant :math:`k`.
    Note that the value of :math:`k` should be chosen with care,
    if its too small, the particle wont follow the biasing potential
    center, if its too large, there will be no overlapping between the 
    different windows.

    The centre of the harmonic potential :math:`x_\text{des}`
    successively takes values from -25 to 25. For each value of
    :math:`x_\text{des}`, an equilibration step of 0.4 ns is
    performed, followed by a step of 2 ns during which the
    position along :math:`x` of the particle is saved in data
    files (one data file per value of :math:`x_\text{des}`). You
    can always increase the duration of the runs for better samplings.

On the choice of k
------------------

..  container:: justify

    As already stated, the difficult part is to choose the value of :math:`k`. You
    want the biasing potential to be strong enough to force
    the atom to move along the axis, and you also want the
    fluctuations of the atom position to be large enough to
    have some overlap in the density probability of two
    neighbor positions, like we have here:

.. figure:: ../figures/level3/free-energy-calculation/overlap-light.png
    :alt: Averaged density profile
    :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/overlap-dark.png
    :alt: Averaged density profile
    :class: only-dark

    Density probability for each run with :math:`k = 1.5` Kcal/mol/Å^2. Note the good
    overlapping between neighbor distributions.

..  container:: justify

    If :math:`k` is too small, the particle never explore the 
    region of interest:

.. figure:: ../figures/level3/free-energy-calculation/overlap015-light.png
    :alt: Averaged density profile
    :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/overlap015-dark.png
    :alt: Averaged density profile
    :class: only-dark

    Density probability for each run with :math:`k = 0.15` Kcal/mol/Å^2. 
    The particle doesn't even explore the middle region.

..  container:: justify

    If :math:`k` is too large, the biasing potential is too large 
    compared to the thermal energy:

.. figure:: ../figures/level3/free-energy-calculation/overlap15-light.png
    :alt: Averaged density profile
    :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/overlap15-dark.png
    :alt: Averaged density profile
    :class: only-dark

    Density probability for each run with :math:`k = 15` Kcal/mol/Å^2. 
    Note the bad overlap between neighbor windows.

WHAM algorithm
--------------

..  container:: justify

    In order to generate the free energy profile from the density distribution, we are going to use
    the WHAM algorithm. You can download and compile the version of |Grossfield|.
    It can be compiled by simply running:

.. |Grossfield| raw:: html

   <a href="http://membrane.urmc.rochester.edu/?page_id=126" target="_blank">Alan Grossfield</a>

..  code-block:: bash

    cd wham
    make clean
    make

..  container:: justify

    It creates an executable called wham that you can 
    copy in the *BiasedSampling/* folder.

    In order to apply the WHAM algorithm to our simulation, we
    first need to create a metadata file. This file simply
    contains the paths of the data files, the value of
    :math:`x_\text{des}`, and the values of :math:`k`. To generate
    the file more easily, you can run this script using Octave
    or Matlab (assuming that the wham algorithm is located in
    the same folder as the LAMMPS simulations):

..  code-block:: bash

    file=fopen('metadata.dat','wt');
    for a=1:50
        X=['./data-k1.5/position.',num2str(a),'.dat ',num2str(a-25),' 1.5'];
        fprintf(file,X);
        fprintf(file,'\n');
    end

..  container:: justify

    The generated file named *metadata.dat* looks like that:

..  code-block:: bash

    ./data-k1.5/position.1.dat -24 1.5
    ./data-k1.5/position.2.dat -23 1.5
    ./data-k1.5/position.3.dat -22 1.5
    ./data-k1.5/position.4.dat -21 1.5
    (...)
    ./data-k1.5/position.48.dat 23 1.5
    ./data-k1.5/position.49.dat 24 1.5
    ./data-k1.5/position.50.dat 25 1.5

..  container:: justify

    Here you can download the |download_metadata| file.
    Then, simply run the following command in the terminal:

.. |download_metadata| raw:: html

   <a href="../../../../../inputs/level3/free-energy-calculation/BiasedSampling/metadata.dat" target="_blank">metadata.dat</a>

..  code-block:: bash

    ./wham -25 25 50 1e-8 119.8 0 metadata.dat PMF.dat

..  container:: justify

    where -25 and 25 are the boundaries, 50 the number of bins,
    1e-8 the tolerance, and 119.8 the temperature. A file named
    PMF.dat has been created, and contains the free energy
    profile in Kcal/mol.

**Results**

..  container:: justify

    We can compare the PMF with we the imposed potential, and
    the agreement in again quite good despite the very short
    calculation time:

.. figure:: ../figures/level3/free-energy-calculation/freeenergy-light.png
    :alt: Result of the umbrella sampling
    :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/freeenergy-dark.png
    :alt: Result of the umbrella sampling
    :class: only-dark

    Calculated potential using umbrella sampling compared to the imposed potential.
    The calculated potential is in blue.

.. include:: ../../contact/accessfile.rst

Going further with exercises
============================

.. include::  ../../contact/requestsolution.rst
        
Monte Carlo versus molecular dynamics
-------------------------------------

..  container:: justify

    Use a Monte Carlo procedure to equilibrate the system
    instead of molecular dynamics. 
    
    Is it more computationally efficient than molecular dynamics?

.. admonition:: Hints (click to reveal)
    :class: dropdown

    Monte Carlo displacement can be made using the fix gcmc command.

Binary fluid
------------

..  container:: justify

    Create a molecular simulation with two species, and apply a
    different potential on them using addforce to create the
    following situation, where one species is mainly trapped in 
    the central part, while the other is mainly trapped in the 
    external part:

.. figure:: ../figures/level3/free-energy-calculation/exercice2-light.png
    :alt: Particles separated by a potential
    :class: only-light

.. figure:: ../figures/level3/free-energy-calculation/exercice2-dark.png
    :alt: Particles separated by a potential
    :class: only-dark

..  container:: justify

    Optional: use particle swapping to exchange between the 2 populations and reach 
    equilibrium faster than with simple molecular dynamics.

Adsorption of ethanol at a solid wall
-------------------------------------

..  container:: justify

    Apply umbrella sampling to calculate the free energy profile
    of a molecule of your choice normal to a solid wall.

.. admonition:: Hints (click to reveal)
    :class: dropdown

    Countless molecules can be downloaded online from the Automated Topology Builder (ATB).
    You can make your life simpler by choosing one of the smallest molecule (something like
    CO2, a small alcohol, or even water).

    Any solid wall would do. You can find graphene-based solid or NaCl crystal surface in
    my |Github_repository|.

.. |Github_repository| raw:: html

   <a href="https://github.com/simongravelle/lammps-input-files" target="_blank">Github repository</a>