.. _contact-before-you-start-label:

Before you start
****************

..  container:: justify

    *LAMMPS tutorials* is made of seven tutorials that are
    ordered by increasing difficulty. The first tutorial, :ref:`lennard-jones-label`,
    is meant for LAMMPS absolute beginners. The complexity
    of the simulations is then progressively increased
    for :ref:`carbon-nanotube-label`,
    :ref:`all-atoms-label`,
    :ref:`sheared-confined-label`,
    and :ref:`reactive-silicon-dioxide-label`.
    Finally, in :ref:`gcmc-silica-label` and
    :ref:`umbrella-sampling-label`,
    some more advanced simulation methods
    are used, namely grand canonical Monte Carlo simulations and
    a free energy method named umbrella sampling.

Required software
=================

LAMMPS (2Aug2023)
-------------------

..  container:: justify

    Download and install LAMMPS version 2Aug2023 by following the
    instructions of the |LAMMPS website installation| :cite:`thompson2022lammps`.
    Depending on your operating system (i.e. Linux, macOS, or Windows),
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

VMD (optional)
--------------

..  container:: justify

    To visualize the simulation, |VMD| version 1.9.3 will
    be used :cite:`humphrey1996vmd`. Some basic instructions for VMD are given here in the
    :ref:`vmd-label`. Feel free to use an alternative visualization
    software like |Ovito|.
    
.. |VMD| raw:: html

   <a href="https://www.ks.uiuc.edu/Research/vmd" target="_blank">VMD</a>
    
.. |Ovito| raw:: html

   <a href="https://www.ovito.org" target="_blank">Ovito</a>
    
Python (optional)
-----------------

..  container:: justify

    To perform post-mortem analysis of the data during the :ref:`mda-label`,
    MDAnalysis version 2.6.1 is used
    together with Python version 3.11.4
    :cite:`van1995python, michaud2011mdanalysis, gowers2016mdanalysis`.

..  container:: justify

    To plot the results from the simulations,
    |Matplotlib Pyplot| version 3.5.2 is used
    in combination with |lammps_logfile|, a library allowing
    one to read the *log* file produced by LAMMPS :cite:`hunter2007Matplotlib, sveinsson2021logfile`.

.. |Matplotlib Pyplot| raw:: html

   <a href="https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html" target="_blank">Matplotlib Pyplot</a>

.. |lammps_logfile| raw:: html

   <a href="https://github.com/henriasv/lammps-logfile" target="_blank">lammps logfile</a>

Text editing software
---------------------

..  container:: justify

    To write and edit LAMMPS input files, a text editor is required.
    Any text editor will do, such as |gedit|, |vim|,
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

Recommended reading
===================

..  container:: justify

    To better understand molecular dynamics simulations, I recommend the reading
    of *Understanding molecular simulation* by Daan Frenkel and Berend
    Smit :cite:`frenkel2023understanding`, as well as
    *Computer simulation of liquids* by Michael Allen and Dominic Tildesley
    :cite:`allen2017computer`. To understand the basic concepts 
    of fluid and Soft Matter systems, I recommend reading *Basic concepts for
    simple and complex liquids* by Jean-Louis Barrat and Jean-Pierre Hansen
    :cite:`barrat2003basic`,
    as well as *Theory of simple liquids: with applications to soft matter*
    by Jean-Pierre Hansen and Ian Ranald McDonald :cite:`hansen2013theory`.