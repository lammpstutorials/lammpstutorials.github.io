.. _contact-before-you-start-label:

Before you start
****************

..  container:: justify

    *LAMMPS tutorials* is made of seven tutorials that are
    ordered by increasing difficulty. :ref:`lennard-jones-label`
    is meant for absolute LAMMPS and molecular dynamics
    beginners. The complexity of the molecular simulation
    is then progressively increased
    for :ref:`carbon-nanotube-label`, 
    :ref:`all-atoms-label`,
    :ref:`sheared-confined-label`,
    and :ref:`reactive-silicon-dioxide-label`.
    Finally, in :ref:`gcmc-silica-label` and
    :ref:`umbrella-sampling-label`,
    some more advanced simulation methods
    are used, namely grand canonical Monte Carlo simulations and 
    a free energy method named umbrella sampling.

Required softwares
==================

..  container:: justify

    The 2 Aug 2023 version of LAMMPS is required
    to follow the tutorials :cite:`thompson2022lammps`.
    The other softwares listed here are optional.

LAMMPS (2 Aug 2023)
-------------------

..  container:: justify

    Download and install the 2 Aug 2023 version of LAMMPS by following the
    instructions of the |LAMMPS website installation|.
    Depending on your operative system (i.e. Linux, macOS, or Windows),
    the procedure may differ.

.. |LAMMPS website installation| raw:: html

   <a href="https://docs.lammps.org/Install.html" target="_blank">LAMMPS website</a>

..  container:: justify

    LAMMPS must be compiled with the following packages:

..  container:: justify

    - MANYBODY
    - MOLECULE
    - KSPACE
    - RIGID
    - REAXFF
    - EXTRA-DUMP

..  container:: justify

    If you decide to use another LAMMPS version, certain commands
    may not work and LAMMPS will throw an |LAMMPS error|.

.. |LAMMPS error| raw:: html

   <a href="https://docs.lammps.org/Errors_messages.html" target="_blank">error message</a>

VMD (1.9.3)
-----------

..  container:: justify

    In order to visualize the simulation, the version 1.9.3 of |VMD| will be used :cite:`humphrey1996vmd`.
    Some basic instructions for VMD are given here in the
    :ref:`vmd-label`. If you prefer, feel free to use an alternative visualization
    software like |Ovito|.
    
.. |VMD| raw:: html

   <a href="https://www.ks.uiuc.edu/Research/vmd" target="_blank">VMD</a>
    
.. |Ovito| raw:: html

   <a href="https://www.ovito.org" target="_blank">Ovito</a>
    
Python (3.11.4)
---------------

..  container:: justify

    To perform post-mortem analysis of the data during the :ref:`mda-label`,
    the version 2.6.1 of MDAnalysis is used
    together with the version 3.11.4 of Python.

..  container:: justify

    In order to plot the results from the simulations,
    the version 3.5.2 of |Matplotlib Pyplot| is used
    in combination with |lammps_logfile|, a library allowing
    one to read the *log* file produced by LAMMPS.

.. |Matplotlib Pyplot| raw:: html

   <a href="https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html" target="_blank">Matplotlib Pyplot</a>

.. |lammps_logfile| raw:: html

   <a href="https://github.com/henriasv/lammps-logfile" target="_blank">lammps logfile</a>

Text editing software
---------------------

..  container:: justify

    In order to write and edit LAMMPS input files, a text editor is required.
    Any text editor will do, such as |gedit|,
    |vim|,
    or |vscode|.
    
.. |gedit| raw:: html

   <a href="https://help.gnome.org/users/gedit/stable/" target="_blank">gedit</a>
    
.. |vim| raw:: html

   <a href="https://www.vim.org/" target="_blank">vim</a>
    
.. |vscode| raw:: html

   <a href="https://code.visualstudio.com/" target="_blank">vscode</a>
    
Find the input scripts
======================

.. include:: ../non-tutorials/accessfile.rst
