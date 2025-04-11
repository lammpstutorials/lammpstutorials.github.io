.. _command-line-label:

Command line
============

LAMMPS can also be executed from the command-line on Linux, macOS, and
Windows without using the GUI.  This is the more common way to run LAMMPS.
Both, the LAMMPS--GUI program and the LAMMPS command-line executable
utilize the same LAMMPS library and thus no changes to the input file are required.

First, open a terminal or command-line prompt window and navigate to the
directory containing the **input.lmp** file.  Then execute:

.. code-block:: bash

    lmp -in input.lmp

where ``lmp`` is the command-line LAMMPS command.

For parallel execution with 4 processors (via OpenMP threads where supported
by the OPENMP package), use:

.. code-block:: bash

    lmp -in input.lmp -pk omp 4 -sf omp

.. admonition:: Note
    :class: non-title-info

    Running in parallel via MPI requires a specially compiled LAMMPS
    package and is not supported by the GUI.  On supercomputers or HPC
    clusters, pre-compiled LAMMPS executables are typically provided
    by the facility's user support team.  For more information, please
    refer to the facility's documentation or contact its user support staff.

See the |LAMMPSdocumentation| for a complete description on how to
run LAMMPS.

.. |LAMMPSdocumentation| raw:: html

    <a href="https://docs.lammps.org/stable/Run_head.html" target="_blank">LAMMPS documentation</a>

