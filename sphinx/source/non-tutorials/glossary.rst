.. _glossary-label:

LAMMPS Glossary
***************

.. container:: justify

    This glossary contains short definitions of the most common LAMMPS-specific words
    encountered when working with LAMMPS. 

Compute
=======

.. container:: justify

    A *compute* is a command that calculates specific quantities during a
    simulation, such as energy, temperature, or pressure. The output(s) of a
    |compute-documentation| named *cname* can be accessed using :math:`\text{c}\_\text{cname}`. 

.. |compute-documentation| raw:: html

   <a href="https://docs.lammps.org/compute.html" target="_blank">computes</a>
    
.. container:: justify

    **Example:** Here, the first compute named *ke_per_atom* returns a vector
    containing the kinetic energy *ke* for each atom, and the second
    compute named *mean_ke* returns the sum of the elements 
    of the vector. The sum is then printed in the log file using *thermo_style*,
    alongside the simulation *step*:

..  code-block:: lammps

    compute ke_per_atom all ke/atom
    compute mean_ke all reduce sum c_ke_per_atom
    thermo_style custom step c_mean_ke

Dump file
=========

.. container:: justify
    
    The *dump* command writes particle coordinates, velocities, and other
    information to a file for post-mortem analysis or visualization.

Ensemble
========

.. container:: justify

    An ensemble refers to a specific thermodynamic state or a collection
    of states that characterize a system of particles. Different ensembles
    represent systems with different degrees of separation from the surrounding
    environment, ranging from completely isolated systems (e.g. microcanonical ensemble)
    to open systems (e.g. grand canonical ensemble).

.. container:: justify

    As an important side note, it is common for LAMMPS beginners to assume that
    their system is in the NVE ensemble (or NVT, or NPT), simply because the
    *fix nve* (or *fix nvt*, or *fix npt*) is used, but this is not necessarily
    the case. For instance, in :ref:`lennard-jones-label`, the input
    contains the following commands:

..  code-block:: lammps

    fix mynve all nve
    fix mylgv all langevin 1.0 1.0 0.1 1530917

.. container:: justify

    Here, the *fix nve* is used to perform the time integration and update
    the position and velocity of the atoms every timestep, and the
    *fix langevin* applies the Langevin thermostat :cite:`schneider1978molecular`. 
    Due to the Langevin thermostat, the system does exchange energy with
    a background implicit solvent, and is therefore not in the NVE (microcanonical)
    ensemble, but rather in the NVT (canonical) ensemble with a
    constant number of particles (N), constant volume (V), and a 
    temperature fluctuating around an equilibrium value (T).

Input script
============

.. container:: justify

    The *input script* is a text file where users define simulation parameters,
    forces, constraints, and other settings using LAMMPS commands.

Fix
===

.. container:: justify

    A *fix* is a command that performs specific tasks during a simulation,
    such as imposing constraints, applying forces, or modifying particle properties.

Force field
===========

.. container:: justify

    A force field is a mathematical model describing the potential energy and forces
    between atoms or particles in a molecular system. It includes the potential energy
    functions (e.g. Lennard-Jones potential) and the parameters (e.g. the 
    Lennard-Jones coefficients :math:`\sigma_{ij}` 
    and :math:`\epsilon_{ij}`).

Grand canonical ensemble
========================

.. container:: justify

    In the grand canonical ensemble, which is also called the macrocanonical ensemble,
    the chemical potential (:math:`\mu`), volume (V), and temperature (T) of the
    system are kept constant. This ensemble is relevant for systems exchanging particles
    with the environment.

LAMMPS
======

.. container:: justify

    LAMMPS is the acronym for Large-scale Atomic/Molecular Massively Parallel Simulator, it
    designs the open-source molecular simulation software.

Minimization
============

.. container:: justify

    Energy minimization refers to the computational process of adjusting the
    atomic positions within a given system to reduce the forces on the atoms
    until they become negligible. Several minimization methods are implemented in LAMMPS,
    including the conjugate gradient :cite:`hestenes1952methods`
    and the steepest descent :cite:`cauchy1847methode`,
    see the |min_style_2| page for an exhaustive list.

.. |min_style_2| raw:: html

    <a href="https://docs.lammps.org/min_style.html" target="_blank">min style</a>

Neighbor list
=============

.. container:: justify

    *Neighbor list* enumerates all pairs of atoms with a separation smaller than
    the cutoff distance.

NVE ensemble
============

.. container:: justify

    In the NVE ensemble, which is also called the microcanonical ensemble,
    the number of particles (N), the volume (V), and the total energy (E) of the
    system are kept constant. This ensemble is relevant for systems that are fully isolated
    and experience no exchange of particles or heat in the environment.

NVT ensemble
============

.. container:: justify

    In the NVT ensemble, which is also called the canonical ensemble,
    the number of particles (N), volume (V), and temperature (T) of the
    system are kept constant. This ensemble is relevant for systems in thermal
    contact with the environment.

NPT ensemble
============

.. container:: justify

    In the NPT ensemble, which is also called the isothermal-isobaric ensemble,
    the number of particles (N), pressure (P), and temperature (T) are kept constant.
    This ensemble is relevant for systems in thermal and mechanical equilibrium
    with the environment.

Pair style
==========

.. container:: justify

    The *pair_style* command sets the potential interactions between pairs of atoms
    (e.g. Lennard-Jones, Coulomb, Morse).

Parallel
========

.. container:: justify

    LAMMPS is designed for parallel computing, enabling simulations to be
    distributed across multiple processors or cores. Parallel execution in
    LAMMPS can be achieved using message-passing parallelism (MPI) or threads.

Restart file
============

.. container:: justify

    A *restart file* allows users to continue simulations from a saved state.

Run
===

.. container:: justify

    A *run* is a command that executes the simulation for a specified number of time steps.

Seed
====

.. container:: justify

    *Seeds* are used to initialize random number generators.
    By setting the seed to a specific value, the user can ensure that the
    sequence of random numbers generated by the simulation will be
    reproducible, which is essential during debugging.

Serial
======

.. container:: justify

    LAMMPS can be run in serial mode on a single processor or core. This is suitable for
    small-scale simulations or when parallel computing resources are not available.

Time step
=========

.. container:: justify

    The simulation progresses through discrete *time steps*.