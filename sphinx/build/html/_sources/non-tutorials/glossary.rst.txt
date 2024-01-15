.. _glossary-label:

Glossary
********

.. container:: justify

    This glossary contains definition of terms commonly encountered when working with LAMMPS. 

Compute
=======

.. container:: justify

    A *compute* is a command that calculates and outputs specific quantities during a
    simulation, such as energy, temperature, or pressure. 

Dump file
=========

.. container:: justify
    
    The *dump* command writes particle coordinates, velocities, and other
    information to a file for post-mortem analysis or visualization.

Ensemble
========

.. container:: justify

    An ensemble refers to a specific thermodynamic state or a collection
    of states that characterizes a system of particles. Different ensembles
    represent systems with different degrees of separation from the surrounding
    environment, ranging from completely isolated systems (e.g. microcanonical ensemble)
    to open systems (e.g. grand canonical ensemble). The ensemble you choose depends on the specific
    problem you are trying to solve and the conditions under which you want to simulate the system.

.. container:: justify

    Careful: It is not because a simulation involves the *fix nve* that the system
    is in the microcanonical or (NVE) ensemble. For instance, in :ref:`lennard-jones-label`,
    the following commands are used:

..  code-block:: lammps

    fix mynve all nve
    fix mylgv all langevin 1.0 1.0 0.1 1530917

.. container:: justify

    Therefore the temperature is imposed by the *fix langevin*,
    and the simulated system is effectively in the NVT ensemble
    (constant number of particles (N), volume (V), and temperature (T)).

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

LAMMPS
======

.. container:: justify

    LAMMPS is the acronym for Large-scale Atomic/Molecular Massively Parallel Simulator, it
    designs the open-source molecular simulation software.

Minimization
============

.. container:: justify

    Energy minimization are used to find a stable configuration for a given system.

Neighbor list
=============

.. container:: justify

    *Neighbor list* enumerates all pairs of atoms with separation smaller than the
    cutoff distance.

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

Serial
======

.. container:: justify

    LAMMPS can be run in serial mode on a single processor or core. This is suitable for
    small-scale simulations or when parallel computing resources are not available.

Time step
=========

.. container:: justify

    The simulation progresses through discrete *time steps*.