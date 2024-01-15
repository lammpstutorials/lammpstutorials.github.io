.. _glossary-label:

Glossary
********

This glossary contains 

Compute
=======

A *compute* is a command that calculates and outputs specific quantities during a
simulation, such as energy, temperature, or pressure. 

Dump file
=========

The *dump* command writes particle coordinates, velocities, and other
information to a file for post-mortem analysis or visualization.

Input script
============

The *input script* is a text file where users define simulation parameters,
forces, constraints, and other settings using LAMMPS commands.

Fix
===

A *fix* is a command that performs specific tasks during a simulation,
such as imposing constraints, applying forces, or modifying particle properties.

Minimization
============

Energy minimization are used to find a stable configuration for a given system.

Neighbor list
=============

*Neighbor list* enumerates all pairs of atoms with separation smaller than the
cutoff distance.

Pair style
==========

The *pair_style* command sets the potential interactions between pairs of atoms
(e.g. Lennard-Jones, Coulomb, Morse).

Parallel
========

LAMMPS is designed for parallel computing, enabling simulations to be
distributed across multiple processors or cores. Parallel execution in
LAMMPS can be achieved using message-passing parallelism (MPI) or threads.

Restart file
============

A *restart file* allows users to continue simulations from a saved state.

Run
===

A *run* is a command that executes the simulation for a specified number of time steps.

Serial
======

LAMMPS can be run in serial mode on a single processor or core. This is suitable for
small-scale simulations or when parallel computing resources are not available.

Time step
=========

The simulation progresses through discrete *time steps*.