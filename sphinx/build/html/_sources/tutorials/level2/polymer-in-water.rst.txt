.. _all-atoms-label:

Polymer in water
****************

.. container:: hatnote

   Solvating and stretching a small polymer molecule

.. figure:: ../figures/level2/polymer-in-water/PEG-dark.webp
    :alt: Movie of a peg molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/PEG-light.webp
    :alt: Movie of a peg molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-light

..  container:: justify

   The goal of this tutorial is to use LAMMPS and solvate a small
   hydrophilic polymer (PEG - PolyEthylene Glycol) in a reservoir of water. 
   An all-atom description is used for all species, and the long
   range Coulomb interaction are solved. Once the system is properly
   equilibrated, a constant stretching force is applied to both
   ends of the polymer, and the evolution of its length with time
   is measured.

..  container:: justify

   This tutorial was inspired by a |Liese2017| by Liese and coworkers, in which
   molecular dynamics simulations are compared with force spectroscopy experiments.

.. |Liese2017| raw:: html

   <a href="https://doi.org/10.1021/acsnano.6b07071" target="_blank">publication</a>

.. include:: ../../non-tutorials/recommand-lj.rst

.. include:: ../../non-tutorials/needhelp.rst

.. include:: ../../non-tutorials/2Aug2023.rst

Preparing water and PEG separately
==================================

..  container:: justify

    In this tutorial, the water is being prepared separately
    from the PEG molecule. The PEG and water will be merged later.

The water
---------

..  container:: justify

    As a first step, a rectangular box of water is created and
    equilibrated at ambient temperature and ambient pressure.
    Create a folder named *pureH2O/*. Inside this folder, create
    an empty text file named *input.lammps*. Copy the following
    lines in it:

..  code-block:: lammps

    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style harmonic
    pair_style lj/cut/coul/long 12
    kspace_style pppm 1e-5
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes

..  container:: justify

    With the unit style *real*, masses are in grams per
    mole, distances in Ångstroms, time in femtoseconds, energies
    in Kcal/mole. With the *atom_style full*, each atom is a dot
    with a mass and a charge that can be
    linked by bonds, angles, dihedrals and impropers. The *bond_style*,
    *angle_style*, and *dihedral_style* commands define the
    potentials for the bonds, angles, and dihedrals used in the simulation,
    here *harmonic*.

..  container:: justify

    Always refer to the LAMMPS |lammps_documentation| if you have doubts about the
    potential used by LAMMPS. For instance, this |lammps_documentation_angle_harmonic|
    gives the expression for the harmonic angular potential.

.. |lammps_documentation| raw:: html

   <a href="https://docs.lammps.org/" target="_blank">documentation</a>

.. |lammps_documentation_angle_harmonic| raw:: html

   <a href="https://docs.lammps.org/angle_harmonic.html" target="_blank">page</a>

..  container:: justify

    Finally, the *special_bonds* command cancels the
    Lennard-Jones interactions between the closest
    atoms of a same molecule.

.. admonition:: About *special bonds*
    :class: info

    Usually, molecular dynamics force fields are parametrized assuming that
    the first neighbors within a molecule do not
    interact directly though LJ or Coulomb potential. Here, since we
    use *lj 0.0 0.0 0.5* and *coul 0.0 0.0 1.0*, the first and second
    neighbors in a molecule only interact through direct bond interactions.
    For the third neighbor (here third neighbor only concerns the PEG molecule, not the water),
    only half of the LJ interaction will be taken into account, and the full Coulomb interaction
    will be used.   

..  container:: justify

    With the *pair_style* named *lj/cut/coul/long*, atoms
    interact through both a Lennard-Jones (LJ) potential and
    Coulombic interactions. The value of :math:`12\,\text{Å}` is 
    the cutoff.

.. admonition:: About cutoff in molecular dynamics
    :class: info

    The cutoff of :math:`12\,\text{Å}` applies to both LJ and Coulombic
    interactions, but in a different way. For LJ *cut*
    interactions, atoms interact with each others only if they
    are separated by a distance smaller than the cutoff. For
    Coulombic *long*, interactions between atoms closer than
    the cutoff are computed directly, and interactions between
    atoms outside that cutoff are computed in the reciprocal space.

..  container:: justify

    Finally the kspace command defines the long-range solver for the (long)
    Coulombic interactions. The *pppm* style refers to
    particle-particle particle-mesh.

.. admonition:: About PPPM
    :class: info

    From |Luty and van Gunsteren|: *'The PPPM method is based on separating the total interaction
    between particles into the sum of short-range
    interactions, which are computed by direct
    particle-particle summation, and long-range interactions,
    which are calculated by solving Poisson's equation using
    periodic boundary conditions (PBCs).'* 

.. |Luty and van Gunsteren| raw:: html

   <a href="https://doi.org/10.1021/jp9518623" target="_blank">Luty and van Gunsteren</a>

..  container:: justify

    Then, let us create a 3D simulation box of dimensions :math:`3 \times 3 \times 3 \; \text{nm}^3`,
    and make space for 8 atom types (1 and 2 for
    the water molecule, and 3 to 8 for the PEG molecule), 6 bond types, 9
    angle types, and 14 dihedrals types.
    Copy the following lines into *input.lammps*:

..  code-block:: lammps

    region box block -15 15 -15 15 -15 15
    create_box 8 box &
    bond/types 6 &
    angle/types 9 &
    dihedral/types 3 &
    extra/bond/per/atom 2 &
    extra/angle/per/atom 1 &
    extra/special/per/atom 2

.. admonition:: About extra per atom commands
    :class: info

    The *extra/x/per/atom* commands are here for
    memory allocation. These commands ensure that enough memory space is left for a
    certain number of attribute for each atom. We wont worry
    about those commands in this tutorial, just keep that in mind if one day you see the following
    error message *ERROR: Molecule topology/atom exceeds system topology/atom*.

..  container:: justify

    Let us create a *PARM.lammps* file containing all the
    parameters (masses, interaction energies, bond equilibrium
    distances, etc). In *input.lammps*, add the following line:

..  code-block:: lammps

    include ../PARM.lammps

..  container:: justify

   Then, download and save the |PARM_PEG.data| file
   next to the *pureH2O/* folder.

.. |PARM_PEG.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/PARM.lammps" target="_blank">parameter</a>

..  container:: justify

    Within *PARM.lammps*, the *mass* and *pair_coeff* of atoms
    of type 1 and 2 are for water, while 
    those of atoms type 3 to 8 are for the PEG
    molecule. Similarly, the *bond_coeff 1* and 
    *angle_coeff 1* are for water, while all
    the other parameters are for the PEG.

..  container:: justify

    Let us create water molecules. To do so, let us
    define what a water molecule is using a molecule *template* called
    *FlexibleH2O.txt*, and then randomly create 350 molecules.
    Add the following lines into *input.lammps*:

..  code-block:: lammps

    molecule h2omol FlexibleH2O.txt
    create_atoms 0 random 350 45615 NULL mol h2omol 14756 overlap 1 maxtry 50

..  container:: justify

    The *overlap 1* option of the *create_atoms* command ensures that no atoms are
    placed exactly at the same position, as this would cause the simulation to
    crash. The *maxtry 50* asks LAMMPS to try at most
    50 times to insert the molecules, which is useful in case some
    insertion attempts are rejected due to overlap. In some case, depending on
    the system and on the values of *overlap*
    and *maxtry*, LAMMPS may not create the desired number of molecules.
    Always check the number of created atoms in the *log* file after
    starting the simulation:

..  code-block:: bw

    Created 1050 atoms

..  container:: justify

    When LAMMPS fails to create the desired number of molecules, a WARNING
    appears in the *log* file.

..  container:: justify

    The molecule template named *FlexibleH2O.txt*
    can be |download_FlexibleH2O.txt|
    and saved in the *pureH2O/* folder.
    This template contains the necessary structural
    information of a water molecule, such as the number of atoms, 
    the id of the atoms that are connected by bonds, by angles, etc.

.. |download_FlexibleH2O.txt| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/pureH2O/FlexibleH2O.txt" target="_blank">downloaded</a>

..  container:: justify

    Then, let us organize the atoms of type 1 and 2 of the water molecules in a group named
    *H2O*, and then perform a small energy minimization. The energy minimization
    is mandatory here given the small *overlap* value of 1 Angstrom chosen in the *create_atoms*
    command. Add the following lines to *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    minimize 1.0e-4 1.0e-6 100 1000
    reset_timestep 0

..  container:: justify

    The *reset_timestep* command is optional. It is used here 
    because the *minimize* command is usually performed over an 
    arbitrary number of steps.

..  container:: justify

    Let us use the *fix npt* to
    control both the temperature and the pressure of the system,
    by adding the following line to *input.lammps*:

..  code-block:: lammps

    fix mynpt all npt temp 300 300 100 iso 1 1 1000

..  container:: justify

    The *fix npt* allows us to impose both a temperature of :math:`300\,\text{K}`
    (with a damping constant of :math:`100\,\text{fs}`),
    and a pressure of 1 atmosphere (with a damping constant of :math:`1000\,\text{fs}`).
    With the *iso* keyword, the three dimensions of the box will be re-scaled simultaneously.

..  container:: justify

    Let us print the atom positions in a *.lammpstrj* file every 1000
    steps (i.e. 1 ps), print the temperature volume, and
    density every 100 steps in 3 separate data files, and
    print the information in the terminal every 1000 steps:

..  code-block:: lammps

   dump mydmp all atom 1000 dump.lammpstrj
   variable mytemp equal temp
   variable myvol equal vol
   fix myat1 all ave/time 10 10 100 v_mytemp file temperature.dat
   fix myat2 all ave/time 10 10 100 v_myvol file volume.dat
   variable myoxy equal count(H2O)/3
   variable mydensity equal ${myoxy}/v_myvol
   fix myat3 all ave/time 10 10 100 v_mydensity file density.dat
   thermo 1000

..  container:: justify

    The variable *myoxy* corresponds to the number of atoms
    divided by 3, i.e. the number of molecules.

.. admonition:: On calling variables in LAMMPS
    :class: info

    Both dollar sign (:math:`\$`)
    and underscore (:math:`v_\text{var}`) can be used to
    call a previously defined variable. 
    With the dollar sign, the initial value of of the variable is returned,
    while with the underscore, the instantaneous value of the variable is returned. 
    To probe the temporal evolution of a variable with time,
    the underscore must be used.

..  container:: justify

    Finally, let us set the timestep to 1.0 fs,
    and run the simulation for 50 ps by adding the
    following lines to *input.lammps*:

..  code-block:: lammps

    timestep 1.0
    run 50000

    replicate 3 1 1

    write_data H2O.data

..  container:: justify

    The *replicate* command is used to replicate the final cubic
    box three times along 
    the *x* direction, thus creating a rectangular box of water 
    just before the final state is written into *H2O.data*.

..
    TO BE PLACED IN A DEDICATED ANNEXE
    .. admonition:: Running LAMMPS in parallel
        :class: info

        This simulation may be a bit slow to complete on 1 single CPU core.
        You can speed it up by running LAMMPS on 2, 4 or even more cores, by typing:

        .. code-block:: bw

            mpirun -np 8 lmp -in input.lammps

        Here, 8 CPU cores are used. The command may vary, depending
        on your OS and LAMMPS installation.

    .. admonition:: Choosing the right number of cores
        :class: dropdown

        When running a simulation in parallel using more than one CPU core, LAMMPS divides the system into
        blocks, and each core is assigned to a given block. Here, as can be seen
        from the terminal, when using *mpirun -np 4*, LAMMPS divides 
        the system into 4 blocks along the x axis:

        .. code-block:: bw
        
            Created orthogonal box = (-40.000000 -15.000000 -15.000000) to (40.000000 15.000000 15.000000)
            4 by 1 by 1 MPI processor grid

        You can force LAMMPS to divide the system differently, let us say along both y and z axis,
        by using the command (in the *input.lammps* file):

        ..  code-block:: bw

            processors 1 2 2

        However, communication between the different cores slows down the computation, so ideally you want 
        to minimize the size of the surface between domains. Here the default choice of LAMMPS (i.e. processors 4 1 1)
        is certainly a better choice.

        If you don't know what is the best number of processors or the best way to cut the system, just perform 
        a short simulation and look at the log file. For instance, if I run the simulation on 1 core I get : 

        .. code-block:: bw

            Performance: 31.567 ns/day, 0.760 hours/ns, 182.680 timesteps/s

        On 4 cores (keeping the default processors 4 1 1):

        .. code-block:: bw

            Performance: 109.631 ns/day, 0.219 hours/ns, 634.440 timesteps/s

        This is much faster, but this is not 4 times faster, because of the cost of communicating between processors.

        On 4 cores and enforcing the stupid choice: processors 1 2 2, I get

        .. code-block:: bw

            Performance: 99.864 ns/day, 0.240 hours/ns, 577.919 timesteps/s

        Its not so bad but still not at good as 4 1 1. 

        On 8 cores (the max I got), I get :

        .. code-block:: bw

            4 by 1 by 2 MPI processor grid

            ...

            Performance: 152.106 ns/day, 0.158 hours/ns, 880.243 timesteps/s

        So LAMMPS chooses to divide the system once along z, and 4 times along x, and the speed is improved, but again the improvement 
        is not linear with the number of cores. 

        **Choose carefully the best number of cores for your simulation so that you don't waste computational resource.
        Sometimes it is better to run 2 simulations on 2 cores each than 1 simulation on 4 cores.**

..  container:: justify

   If you open the *dump.lammpstrj* file using VMD, you should
   see the system quickly reaching its equilibrium volume and density.

.. figure:: ../figures/level2/polymer-in-water/water-light.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-light

.. figure:: ../figures/level2/polymer-in-water/water-dark.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-dark

.. container:: figurelegend

    Figure: Water reservoir after equilibration and after 
    the *replicate 3 1 1* command. Oxygen atoms are in red, and 
    hydrogen atoms in white.

..  container:: justify

   You can also open the *density.dat* file
   to ensure that the system converged toward an equilibrated
   liquid water system during the 50 ps of simulation.

.. figure:: ../figures/level2/polymer-in-water/density_H2O-light.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-light

.. figure:: ../figures/level2/polymer-in-water/density_H2O-dark.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-dark

.. container:: figurelegend

    Figure: Evolution of the density of water with time. The
    density :math:`\rho` reaches
    a plateau after :math:`\approx 30\,\text{ps}`.

..  container:: justify

   If needed, you can |download_H2O.data|
   the water reservoir I have equilibrated and continue with
   the tutorial.

.. |download_H2O.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/pureH2O/H2O.data" target="_blank">download</a>

The PEG molecule
----------------

..  container:: justify

    Now that the water box is ready, let us prepare the PEG
    molecule in an empty box. Create a second folder next to *pureH2O/*, call it
    *singlePEG/*, and create a new blank file called *input.lammps*
    in it. Copy the same first lines as previously:

..  code-block:: lammps

    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style harmonic
    pair_style lj/cut/coul/long 12
    kspace_style pppm 1e-5
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes

..  container:: justify

    Let us read the original positions for the atoms of the PEG molecule, as
    well as the same parameter file as previously, by adding the following 
    lines to *input.lammps*:

..  code-block:: lammps

    read_data init.data
    include ../PARM.lammps

..  container:: justify

   |download_init.data|
   the init.data file and save it in the *singlePEG/* folder.
   It contains the initial parameters of the PEG molecules
   (atoms, bonds, charges, etc.) that was downloaded from the ATB repository.

.. |download_init.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/singlePEG/init.data" target="_blank">Download</a>

..  container:: justify

   Let us print the atom positions and thermodynamic
   information very frequently (because we anticipate that the
   energy minimization will be short). Add the following lines 
   to *input.lammps*:

..  code-block:: lammps

    dump mydmp all atom 10 dump.lammpstrj
    thermo 1

..  container:: justify

    Next, let us perform a minimization of energy. Here, this
    step is required because the initial configuration of the
    PEG molecule is really far from equilibrium. Add the following line
    to *input.lammps*:

..  code-block:: lammps

    minimize 1.0e-4 1.0e-6 100 1000

..  container:: justify

    After the minimization, the high frequency dump command is
    cancelled, and a new dump command with lower frequency is
    used (see below). We also reset the time to 0 with
    *reset_timestep* command:

..  code-block:: lammps

    undump mydmp
    reset_timestep 0

..  container:: justify

    The PEG is then equilibrated in the NVT ensemble (fix NVE +
    temperature control = NVT). No box relaxation is required as
    the PEG is in vacuum:

..  code-block:: lammps

    fix mynve all nve
    fix myber all temp/berendsen 300 300 100

..  container:: justify

    Let us print the temperature in a file by also 
    adding the following lines to *input.lammps*:

..  code-block:: lammps

    dump mydmp all atom 1000 dump.lammpstrj
    dump_modify mydmp append yes
    thermo 1000
    variable mytemp equal temp
    fix myat1 all ave/time 10 10 100 v_mytemp file output-temperature.dat

..  container:: justify

    The *dump_modify* ensures that the coordinates are written 
    in the already existing *dump.lammpstrj* file. 
    Finally let us run the simulation for a very short time (10 ps)
    by adding the following lines to *input.lammps*:

..  code-block:: lammps

    timestep 1.0
    run 10000
    write_data PEG.data

..  container:: justify

    If you open the *dump.lammpstrj* file
    using VMD, you can see the PEG molecule 
    gently equilibrating until reaching a reasonable state.

.. figure:: ../figures/level2/polymer-in-water/singlePEG-light.png
    :alt: PEG in vacuum as simulated with LAMMPS
    :class: only-light

.. figure:: ../figures/level2/polymer-in-water/singlePEG-dark.png
    :alt: PEG in vacuum as simulated with LAMMPS
    :class: only-dark

..  container:: figurelegend

    Figure: The PEG molecule in vacuum.
    The carbon atoms are in pink, the oxygen atoms in red, and the hydrogen
    atoms in white. 
    
..  container:: justify

   Alternatively, you can |download_PEG.data|
   the PEG molecule I have equilibrated and continue with the tutorial.

.. |download_PEG.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/singlePEG/PEG.data" target="_blank">download</a>

Solvated PEG
============

..  container:: justify

    Once both the water and the PEG are equilibrated, 
    we can safely merge the two systems before performing the 
    pull experiment on the polymer.

Mixing the PEG with water
-------------------------

..  container:: justify

    To merge the two systems, let us:

    - (1) import both previously generated data files (PEG.data and H2O.data) into the same simulation,
    - (2) delete the overlapping molecules, and 
    - (3) re-equilibrate the new system. 

..  container:: justify

    Create a third folder alongside *pureH2O/* and *singlePEG/*,
    and call it *mergePEGH2O/*. Create a new blank file in it,
    called *input.lammps*. Within *input.lammps*, copy the same first lines as
    previously:

..  code-block:: lammps

    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style harmonic
    pair_style lj/cut/coul/long 12
    kspace_style pppm 1e-5
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes

..  container:: justify

    Then, import the two previously generated data files *H2O.data*
    and *PEG.data*, as well as the *PARM.lammps* file
    by adding that to *input.lammps*:

..  code-block:: lammps

    read_data ../pureH2O/H2O.data extra/bond/per/atom 2 extra/angle/per/atom 5 extra/dihedral/per/atom 9
    read_data ../singlePEG/PEG.data add append shift 25 0 0
    include ../PARM.lammps

..  container:: justify

    When using the *read_data* command more than once, one needs
    to use the *add append* keyword. When doing so, the
    simulation box is initialized by the first *read_data* only, and the 
    second *read_data* only imports additional atoms.

..  container:: justify

    The *extra/x/per/atom* commands are again here for memory allocation.
    The *shift 25 0 0* that is applied to the polymer is there
    to recenter the polymer in the rectangular box by shifting its position 
    by 25 Angstroms along the *x* axis.

..  container:: justify

    Let us create 2 groups to differentiate the PEG from the H2O,
    by adding the following lines to *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    group PEG type 3 4 5 6 7 8

..  container:: justify

    Water molecules that are overlapping with the PEG must be
    deleted to avoid future crashing. Add the following line 
    to *input.lammps*:

..  code-block:: lammps

    delete_atoms overlap 2.0 H2O PEG mol yes

..  container:: justify

    Here, the value of 2 Angstroms for the overlap cutoff was fixed arbitrarily,
    and can be chosen through trial and error. If the cutoff is too small, the 
    simulation will crash. If the cutoff it too large, too many water molecules
    will unnecessarily be deleted.

..  container:: justify

    Finally, let us use the *fix NPT* to control the
    temperature, as well as the pressure by allowing the 
    box size to be rescalled along the *x* axis:

..  code-block:: lammps

    fix mynpt all npt temp 300 300 100 x 1 1 1000
    timestep 1.0

..  container:: justify

    Once more, let us dump the atom positions and a few
    information about the evolution of the simulation:

..  code-block:: lammps

    dump mydmp all atom 100 dump.lammpstrj
    thermo 100
    variable mytemp equal temp
    variable myvol equal vol
    fix myat1 all ave/time 10 10 100 v_mytemp file temperature.dat
    fix myat2 all ave/time 10 10 100 v_myvol file volume.dat

..  container:: justify

    Let us also print the total enthalpy:

..  code-block:: lammps

    variable myenthalpy equal enthalpy
    fix myat3 all ave/time 10 10 100 v_myenthalpy file enthalpy.dat

..  container:: justify

    Finally, let us perform a short equilibration and print the
    final state in a data file. Add the following lines to the 
    data file:

..  code-block:: lammps

    run 10000
    write_data mix.data

..  container:: justify

    If you open the *dump.lammpstrj* file using VMD, 
    or have a look at the evolution of the volume in *volume.dat*,
    you should see that the box dimension slightly evolves along *x*
    to accomodate the new configuration.
    
..  container:: justify

    The final system looks like that:

.. figure:: ../figures/level2/polymer-in-water/solvatedPEG_light.png
   :alt: PEG in water
   :class: only-light

.. figure:: ../figures/level2/polymer-in-water/solvatedPEG_dark.png
   :alt: PEG in water
   :class: only-dark

.. container:: figurelegend

   Figure: A single PEG molecule in water. 
   Some water molecules are represented as a transparent continuum 
   field for clarity.

Stretching the PEG molecule
---------------------------

..  container:: justify

   Here, a constant forcing is applied to the two ends of the
   PEG molecule until it stretches. Create a new folder next
   to the 3 previously created folders, call it *pullonPEG/*
   and create a new input file in it called *input.lammps*.

..  container:: justify

   First, let us create a variable *f0* corresponding to the magnitude
   of the force we are going to apply. The force magnitude is
   chosen to be large enough to overcome the thermal
   agitation and the entropic contribution from both water
   and PEG molecules (it was chosen by trial and error). Copy
   in the *input.lammps* file:

..  code-block:: lammps

    variable f0 equal 5 # kcal/mol/A # 1 kcal/mol/A = 67.2 pN

..  container:: justify

    Then, as previouly, copy:

..  code-block:: lammps

    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style harmonic
    pair_style lj/cut/coul/long 12
    kspace_style pppm 1e-5
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes

..  container:: justify

    Start the simulation from the equilibrated PEG and water
    system, and include again the parameter file by
    adding to the *input.lammps*:

..  code-block:: lammps

    read_data ../mergePEGH2O/mix.data
    include ../PARM.lammps

..  container:: justify

    Then, let us create 4 atom groups: H2O and PEG (as
    previously) as well as 2 groups containing only the 
    2 oxygen atoms of types 6 and 7, respectively.
    Atoms of types 6 and 7 correspond to the oxygen atoms located at the
    ends of the PEG molecule, which we are going to use 
    to pull on the PEG molecule. Add the following lines to
    the *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    group PEG type 3 4 5 6 7 8
    group topull1 type 6
    group topull2 type 7

..  container:: justify

    Let us add the *dump* command again to print the atom positions:

..  code-block:: lammps

   dump mydmp all atom 1000 dump.lammpstrj

..  container:: justify

    Let us use a simple thermostating for all atoms by adding the 
    following lines to *input.lammps*:

..  code-block:: lammps

    timestep 1.0
    fix mynvt all nvt temp 300 300 100

..  container:: justify

    Let us also print the end-to-end distance of the PEG,
    here defined as the distance between the groups *topull1*
    and *topull2*, as well as the temperature of the system 
    by adding the following lines to *input.lammps*:

..  code-block:: lammps

    variable mytemp equal temp
    fix myat1 all ave/time 10 10 100 v_mytemp file output-temperature.dat
    variable x1 equal xcm(topull1,x)
    variable x2 equal xcm(topull2,x)
    variable y1 equal xcm(topull1,y)
    variable y2 equal xcm(topull2,y)
    variable z1 equal xcm(topull1,z)
    variable z2 equal xcm(topull2,z)
    variable delta_r equal sqrt((v_x1-v_x2)^2+(v_y1-v_y2)^2+(v_z1-v_z2)^2)
    fix myat2 all ave/time 10 10 100 v_delta_r file output-end-to-end-distance.dat
    thermo 1000

..  container:: justify

    Finally, let us run the simulation for 30 ps without
    any external forcing:

..  code-block:: lammps

    run 30000

..  container:: justify

    This first run serves a benchmark to quantify the changes
    induced by the forcing. Then, let us apply a forcing on the 2 oxygen
    atoms using two *add_force* commands, and run for an extra 30 ps:

..  code-block:: lammps

    fix myaf1 oxygen_end1 addforce ${f0} 0 0
    fix myaf2 oxygen_end2 addforce -${f0} 0 0
    run 30000

..  container:: justify

    If you open the *dump.lammpstrj* file using *VMD*, you should
    see that the PEG molecule eventually aligns in the direction
    of the force:

.. figure:: ../figures/level2/polymer-in-water/pulled_peg_dark.png
    :alt: PEG molecule in water
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/pulled_peg_light.png
    :alt: PEG molecule in water
    :class: only-light

.. container:: figurelegend

    Figure: PEG molecule streched along the *x* direction in water.
    Some water molecules are represented as a transparent continuum 
    field for clarity.

..  container:: justify

    The evolution of the end-to-end
    distance over time shows the PEG adjusting
    to the external forcing:

.. figure:: ../figures/level2/polymer-in-water/distance-dark.png
    :alt: plot of the end-to-end distance versus time
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/distance-light.png
    :alt: plot of the end-to-end distance versus time
    :class: only-light

.. container:: figurelegend

    Figure: Evolution of the end-to-end distance of the PEG molecule
    with time. The forcing starts at :math:`t = 30` ps.

.. include:: ../../non-tutorials/accessfile.rst

Going further with exercises
============================

.. include:: ../../non-tutorials/link-to-solutions.rst

Extract radial distribution function
------------------------------------

..  container:: justify

    Extract the radial distribution functions (RDF or :math:`g(r)`)
    between water molecules, as well as between water and PEG molecules:

.. figure:: ../figures/level2/polymer-in-water/RDF-dark.png
    :alt: RDF g(r) for water and peg
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/RDF-light.png
    :alt: RDF g(r) for water and peg
    :class: only-light

.. container:: figurelegend

    Figure: Radial distribution function between the oxygen atoms 
    of water, as well as between the oxygen atoms of water and the 
    carbon and oxygen atoms of the PEG molecule.  

Add salt to the mixture
-----------------------

..  container:: justify

    Add some :math:`\text{Na}^+` and 
    :math:`\text{Cl}^-` ions to the mixture using the method
    of your choice. :math:`\text{Na}^+` ions are 
    characterised by their mass :math:`m = 22.98\,\text{g/mol}`,
    their charge :math:`q = +1\,e`, and Lennard-Jones
    parameters, :math:`\epsilon = 0.0469\,\text{kcal/mol}, \sigma = 0.243\,\text{nm}`,
    and :math:`\text{Cl}^-` ions by their
    mass :math:`m = 35.453\,\text{g/mol}`,
    charge :math:`q = -1\,e` and Lennard-Jones
    parameters, :math:`\epsilon = 0.15\,\text{kcal/mol}, \sigma = 0.4045\,\text{nm}`.

.. figure:: ../figures/level2/polymer-in-water/salt-exercise-dark.png
    :alt: PEG in a NaCl solution
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/salt-exercise-light.png
    :alt: PEG in a NaCl solution
    :class: only-light

.. container:: figurelegend

    Figure: A PEG molecule in water, with :math:`\text{Na}^+` ions in 
    purple and :math:`\text{Cl}^-` ions in cyan.
