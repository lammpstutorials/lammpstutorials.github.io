.. _all-atoms-label:

Polymer in water
****************

.. container:: hatnote

   Solvating and stretching a small polymer molecule

.. figure:: ../figures/level2/polymer-in-water/video-PEG-dark.webp
    :alt: Movie of a peg molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-dark

.. figure:: ../figures/level2/polymer-in-water/video-PEG-light.webp
    :alt: Movie of a peg molecule in water as simulated with LAMMPS
    :height: 250
    :align: right
    :class: only-light

..  container:: justify

   The goal of this tutorial is to use LAMMPS and
   create a small hydrophilic polymer (PEG -
   PolyEthylene Glycol) in a reservoir of water. 
   An all-atom description is used for all species.
   Once the system is created, a constant stretching force is applied to both
   ends of the polymer, and the evolution of its length with time
   will be measured.

..  container:: justify

   This tutorial was inspired by a very nice |Liese2017| by Liese and coworkers, in which
   they compare molecular dynamics simulations with force spectroscopy experiments.

.. |Liese2017| raw:: html

   <a href="https://doi.org/10.1021/acsnano.6b07071" target="_blank">publication</a>

.. include:: ../../contact/recommand-lj.rst

.. include:: ../../contact/needhelp.rst

.. include:: ../../contact/2Aug2023.rst

Preparing water and PEG separately
==================================

..  container:: justify

    As for most simulations, several possible routes can be used
    to create the system. In this tutorial, the water is being prepared separately
    from the PEG molecule. Then, PEG and water are merged.

The water
---------

..  container:: justify

    As a first step, a rectangular box of water is created and
    equilibrated at ambient temperature and ambient pressure.
    Create a folder named *pureH2O/*. Inside this folder, create
    an empty text file named *input.lammps*. Copy the following
    lines in it:

..  code-block:: lammps

    # LAMMPS input script
    units real
    atom_style full
    bond_style harmonic
    angle_style harmonic
    dihedral_style harmonic
    pair_style lj/cut/coul/long 12
    kspace_style pppm 1e-5
    special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes

..  container:: justify

    With the unit style *real*,
    masses are in grams per
    mole, distances in Ångstroms, time in femtoseconds, energies
    in Kcal/mole. With the *atom_style full*, each atom is a dot
    with a mass and a charge that can be
    linked by bonds, angles, dihedrals and impropers potentials
    (for example to form molecules). The *bond_style*,
    *angle_style*, and *dihedral_style* commands define the
    styles of bond angle, and dihedrals used in the simulation,
    respectively, and the *harmonic* keyword
    imposes the form of the potentials.

..  container:: justify

    Finally, the *special_bonds* command to cancel the
    Lennard-Jones interactions between the closest
    atoms of a same molecule.

.. admonition:: About *special bonds*
    :class: info

    Usually, molecular dynamics force fields are parametrized assuming that the first neighbors within a molecule do not
    interact directly. Here, since we use *lj 0.0 0.0 0.5*, the first and second neighbors in a molecule won't interact
    with each other through LJ potentials, and therefore they only interact through direct bond interactions.
    For the third neighbor (here third neighbor only concerns the PEG molecule, not the water),
    only half of the LJ interaction will be added.   

..  container:: justify

    With the *pair_style* named *lj/cut/coul/long*, atoms
    interact through both a Lennard-Jones (LJ) potential and
    through Coulombic interactions. The value of *12* is 
    the cut off for the short range interactions.

.. admonition:: About cutoff in molecular dynamics
    :class: info

    The cutoff of 12 Ångstroms applies to both LJ and Coulombic
    interactions, but in a different way. For LJ *cut*
    interactions, atoms interact with each others only if they
    are separated by a distance smaller than the cutoff. For
    Coulombic *long*, interaction between atoms closer than
    the cutoff are computed directly, and interaction between
    atoms outside that cutoff are computed in the reciprocal
    space.

..  container:: justify

    Finally the kspace command defines the long-range solver for the (long)
    Coulombic interactions. The *pppm* style refers to
    particle-particle particle-mesh.

.. admonition:: Background Information (optional) -- About PPPM
    :class: dropdown

    *The PPPM
    method is based on separating the total interaction
    between particles into the sum of short-range
    interactions, which are computed by direct
    particle-particle summation, and long-range interactions,
    which are calculated by solving Poisson's equation using
    periodic boundary conditions (PBCs).* 
    |Luty and van Gunsteren|

.. |Luty and van Gunsteren| raw:: html

   <a href="https://doi.org/10.1021/jp9518623" target="_blank">Luty and van Gunsteren</a>

..  container:: justify

    Then, let us create a 3D simulation box of dimensions :math:`3 \times 3 \times 3 \; \text{nm}^3`,
    and make space for 7 atom types (1 and 2 for
    the water oxygen and hydrogen, respectively, and 3, 4, 5, 6
    and 7 for the PEG molecule (see below)), 6 bond types, 9
    angle types, and 14 dihedrals types.
    Copy the following lines into *input.lammps*:

..  code-block:: lammps

    region box block -15 15 -15 15 -15 15
    create_box 7 box &
    bond/types 6 &
    angle/types 9 &
    dihedral/types 3 &
    extra/bond/per/atom 2 &
    extra/angle/per/atom 1 &
    extra/special/per/atom 2

.. admonition:: About extra per atom commands
    :class: info

    The *extra/something/per/atom* commands are here for
    memory allocation. These commands ensure that enough memory space is left for a
    certain number of attribute for each atom. We wont worry
    about those commands in this tutorial, just keep that in mind if one day you see the following
    error message *ERROR: Molecule topology/atom exceeds system topology/atom*:

..  container:: justify

    Let us create a *PARM.lammps* file containing all the
    parameters (masses, interaction energies, bond equilibrium
    distances, etc). In *input.lammps*, add the following line:

..  code-block:: lammps

    include ../PARM.lammps

..  container:: justify

   Then, next to the *pureH2O/* folder, create a blank file called
   *PARM.lammps* and copy the following lines in it:

..  code-block:: lammps

    mass 1 15.9994 # H2O O
    mass 2 1.008 # H2O H
    mass 3 12.011 # CC32A
    mass 4 15.9994 # OC30A
    mass 5 1.008 # HCA2
    mass 6 15.9994 # OC311
    mass 7 1.008 # HCP1

    pair_coeff 1 1 0.119431704 3.400251
    pair_coeff 2 2 0.0 0.0
    pair_coeff 3 3 0.25265643 2.8491607 
    pair_coeff 4 4 0.06630155 3.5811794 
    pair_coeff 5 5 0.028293679 2.373408  
    pair_coeff 6 6 0.0 0.0 
    pair_coeff 7 7 0.11949714 3.1000042 

    bond_coeff 1 442.1606 0.972 
    bond_coeff 2 1109.2926 1.12 
    bond_coeff 3 399.79163 1.43 
    bond_coeff 4 400.0343 1.53 
    bond_coeff 5 179.2543 0.971 
    bond_coeff 6 155.35373 1.42

    angle_coeff 1 47.555878 103.0 
    angle_coeff 2 30.173132 109.5 
    angle_coeff 3 47.69405 109.5 
    angle_coeff 4 55.113907 111.0 
    angle_coeff 5 65.47197 111.3 
    angle_coeff 6 54.993103 110.3 
    angle_coeff 7 55.0234 111.4 
    angle_coeff 8 180.46019 109.0 
    angle_coeff 9 30.173132 110.0 

    dihedral_coeff 1 0.30114722 1 3 
    dihedral_coeff 2 1.414914 1 3 
    dihedral_coeff 3 0.0 1 1 

..  container:: justify

    The *mass* and *pair_coeff* of atoms of type 1 and 2 are for water, while 
    those of atoms type 3 to 7 are for the PEG. In addition, *bond_coeff 1* and 
    *angle_coeff 1* are for water, while all the other parameters are for the PEG.

..  container:: justify

    Let us create water molecules. To do so, let us
    define what a water molecule is using a molecule *template* called
    *FlexibleH2O.txt*, and then randomly create 350 molecules.
    Add the following lines into *input.lammps*:

..  code-block:: lammps

    molecule h2omol FlexibleH2O.txt
    create_atoms 0 random 350 456415 NULL mol h2omol 454756 overlap 1.0 maxtry 50

..  container:: justify

    The *overlap 1.0* option of the *create_atoms* command ensures that no atom are
    placed exactly at the same position as it would make the simulation
    crashes. The *maxtry 50* asks LAMMPS to try at least
    50 times to insert molecules. In some case, depending on the system and on the *overlap*
    and *maxtry* values, LAMMPS may not create the desired number of molecules. We 
    can check in the *log* file that it is indeed the case:

..  code-block:: bw

    Created 1050 atoms

..  container:: justify

    When LAMMPS fails to create the desired number of molecules, a WARNING can 
    be seen in the *log* file.

..  container:: justify

    The molecule template named *FlexibleH2O.txt*
    can be |download_FlexibleH2O.txt|
    and saved in the *pureH2O/* folder.
    This template contains all the necessary structural
    information of a water molecule, such as the number of atoms, 
    which pair of atoms are connected by bonds, which
    groups of atoms are connected by angles, etc.

.. |download_FlexibleH2O.txt| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/pureH2O/FlexibleH2O.txt" target="_blank">downloaded</a>

..  container:: justify

    Then, let us group the atoms of the water molecules in a group named
    *H2O*, and then perform a small energy minimization because the 
    *overlap* value of 1 Angstrom chosen in the *create_atoms* command 
    is way too small for water, and we may expect the system to fail
    in absence of minimization. Add the following lines to *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    minimize 1.0e-4 1.0e-6 100 1000
    reset_timestep 0

..  container:: justify

    The *reset_timestep* command is optional. It is used here 
    because the *minimize* command is usually performed over an 
    arbitrary number of steps, which is some situation can make
    life difficult in future runs.

..  container:: justify

    Let us use the fix NPT to
    control both the temperature and the pressure,
    by adding the following line to *input.lammps*:

..  code-block:: lammps

    fix mynpt all npt temp 300 300 100 iso 1 1 1000

..  container:: justify

    The fix NPT allows us to impose both a temperature of 300 K (with a damping constant of 100 fs),
    and a pressure of 1 atmosphere (with a damping constant of 1000 fs). With the *iso* keyword, the
    three dimensions of the box will be re-scaled simultaneously.

..  container:: justify

    Let us print the atom positions in a dump file every 1000
    steps (i.e. 1 ps), print the temperature volume, and
    density every 100 steps in 3 separate data files, and
    print the information in the terminal every 1000 steps:

..  code-block:: lammps

   dump mydmp all atom 1000 dump.lammpstrj
   variable mytemp equal temp
   variable myvol equal vol
   fix myat1 all ave/time 10 10 100 v_mytemp file temperature.dat
   fix myat2 all ave/time 10 10 100 v_myvol file volume.dat
   variable myoxy equal count(H2O)/3 # divide by 3 to get the number of molecule, not atom
   variable mydensity equal ${myoxy}/v_myvol
   fix myat3 all ave/time 10 10 100 v_mydensity file density.dat
   thermo 1000

.. admonition:: On calling variables
    :class: info

    Both :math:`\$ \text{var}` and :math:`v_\text{var}` can be used to call a previously defined variable named *var*. 
    However,  :math:`\$ \text{var}` returns the initial value of *var*,
    while :math:`v_\text{var}` returns the instantaneous value of *var*. 

..  container:: justify

    In the formula for the density (number of
    molecule divided by volume), the underscore *_* is used to
    call myvol because the volume is expected
    to evolve in time when using fix NPT, but the dollar sign :math:`\$` is used to call *myoxy* as the
    number of molecules is not expected to evolve during the
    simulation.

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

    The *replicate* command replicate the cubic box three times along 
    the *x* direction, to create a rectangular box of water 
    right before the final state is written into *H2O.data*.

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
   Here is a snapshot of the final state of the simulation, after 
   the *replicate 3 1 1* command:

.. figure:: ../figures/level2/polymer-in-water/water-light.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-light

.. figure:: ../figures/level2/polymer-in-water/water-dark.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-dark

..  container:: justify

   You can also open the *density.dat* file
   to ensure that the system converged toward an equilibrated
   liquid water system during the 50 ps of simulation:

.. figure:: ../figures/level2/polymer-in-water/density_H2O-light.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-light

.. figure:: ../figures/level2/polymer-in-water/density_H2O-dark.png
    :alt: Curves showing the equilibration of the water reservoir
    :class: only-dark

.. container:: figurelegend

    Figure: Evolution of the density of water with time.

..  container:: justify

   If needed, you can |download_H2O.data|
   the water reservoir I have equilibrated and continue with
   the tutorial.

.. |download_H2O.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/pureH2O/H2O.data" target="_blank">download</a>

.. include:: ../../contact/supportme.rst

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

    The PEG molecule in vacuum.
    The carbon atoms are in pink, the oxygen atoms in red, and the hydrogen
    atoms in white. 
    
..  container:: justify

   Alternatively, you can |download_PEG.data|
   the PEG molecule I have equilibrated and continue with the tutorial.

.. |download_PEG.data| raw:: html

   <a href="../../../../../inputs/level2/polymer-in-water/singlePEG/PEG.data" target="_blank">download</a>

Solvating the PEG molecule
==========================

..  container:: justify

    Let us now merge the PEG molecule and the
    water reservoir. We do it by:

    - (1) importing both previously generated data files (PEG.data and H2O.data) into the same simulation,
    - (2) deleting the overlapping molecules, and 
    - (3) re-equilibrating the new system. 

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

    Then, import the two previously generated data files, as well as the
    same parameter file by adding that to *input.lammps*:

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

    The use of the *extra/x/per/atom* commands is only for memory allocation issue,
    and the *shift 25 0 0* applied to the polymer is there to recenter the polymer in
    the rectangular box.

..  container:: justify

    Let us create 2 groups to differentiate the PEG from the H2O,
    by adding the following lines to *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    group PEG type 3 4 5 6 7

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
    temperature, as well as the pressure along the *x* axis:

..  code-block:: lammps

    fix mynpt all npt temp 300 300 100 x 1 1 1000
    timestep 1.0

..  container:: justify

    Once more, let us dump the atom positions and a few
    information about the evolution simulation:

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

   A single PEG molecule in water. 
   Water molecules are represented as cyan sticks for clarity.

Stretching the PEG molecule
===========================

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

    Then, let us create 3 atom groups: H2O and PEG (as
    previously) as well as a groups containing only the 
    2 oxygen atoms of types 6
    corresponding to the oxygen atoms located at the
    ends of the PEG molecule, which we are going to use 
    to pull on the PEG molecule. Add the following lines to
    the *input.lammps*:

..  code-block:: lammps

    group H2O type 1 2
    group PEG type 3 4 5 6 7

    group topull type 6
    write_dump topull atom topull.lammpstrj
    run 0

..  container:: justify

    Here, we use the *write_dump* command to print the 
    *topull* group in a file. Execute the *input.lammps*
    file using LAMMPS so that *topull.lammpstrj* get written.
    In my case, the two last lines of *topull.lammpstrj* reads:

..  code-block:: lammps

    (...)
    3216 6 0.339394 0.0930772 0.737381
    3151 6 0.523881 0.760453 0.520978

..  container:: justify

    From these two lines, one can extract the two ids 
    of the atoms of types 6, 3216 and 3151, respectively. 
    Let us then create two additional groups by adding to 
    *input.lammps*:

..  code-block:: lammps

    group topull1 id 3216
    group topull2 id 3151

..  container:: justify

    The values of the atom id may be different in your case. 
    Be sure to use the id values corresponding to the atoms 
    of your system.

..  container:: justify

    Let us add the *dump* command again to print the atom positions:

..  code-block:: lammps

   dump mydmp all atom 1000 dump.lammpstrj
   # write_dump all atom dump.lammpstrj
   # dump myxtc xtc atom 1000 dump.xtc

.. admonition:: Use less disk space by using the xtc format
    :class: info
    
    To generate smaller dump files, use the
    compressed *xtc* format. You can do it by commenting the
    mydmp line and by uncommenting both the *write_dump* and
    *myxtc* lines. *xtc* files are compressed, and not readable
    by humans, contrarily to the LAMMPS native format *lammpstrj*. 

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

    Finally, let us run the simulation for 10 ps without
    any external forcing:

..  code-block:: lammps

    run 10000

..  container:: justify

    This 10 ps serves as an extra small equilibration. In principle, 
    it is not necessary as equilibration was properly performed during the 
    previous step. Then, let us apply a forcing on the 2 oxygen atoms using 2
    *add_force* commands, and run for an extra 50 ps:

..  code-block:: lammps

    fix myaf1 oxygen_end1 addforce ${f0} 0 0
    fix myaf2 oxygen_end2 addforce -${f0} 0 0
    run 50000

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

    PEG molecule streched along the *x* direction in water.
    Water molecules are represented as cyan sticks for clarity.

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

    Evolution of the end-to-end distance of the PEG molecule
    with time. The forcing starts at :math:`t = 10` ps.

.. include:: ../../contact/accessfile.rst

Going further with exercises
============================

.. include:: ../../contact/requestsolution.rst

Generate a PEG-H2O mixture
--------------------------

..  container:: justify

    Create a PEG-H2O mixture with several PEG molecules hydrated in a
    cubic box.

.. admonition:: Hints
    :class: info

    Have a look at the LAMMPS *replicate* command.
    Note tthat there is no obligation to equilibrate the water molecules separately from the PEG,
    as we did here. You can also create the water molecules directly around the PEG molcule
    using the *create_atom* command.

Post-mortem analysis
--------------------

..  container:: justify

    In today research, most data analyses are
    done after the simulation is over, and it is important for
    LAMMPS users to know how to do it.

    Import the trajectory using Python, and re-extract the
    end-to-end distance.

.. admonition:: Hints
    :class: dropdown

    You can import *lammpstrj* file using *MDAnalysis* in *Python*:

    ..  code-block:: bw

        u = mda.Universe("dump.lammpstrj", format = "LAMMPSDUMP")