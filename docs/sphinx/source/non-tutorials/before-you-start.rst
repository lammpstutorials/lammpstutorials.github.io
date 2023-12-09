.. _contact-before-you-start-label:

Before you start
****************

..  container:: justify

    *LAMMPS tutorials* is made of seven tutorials that are
    ordered by increasing difficulty. :ref:`lennard-jones-label`
    is meant for absolute LAMMPS and molecular dynamics
    beginners, and the complexity of the simulation
    is progressively increased
    for :ref:`carbon-nanotube-label`, 
    :ref:`all-atoms-label`,
    :ref:`sheared-confined-label`,
    and :ref:`reactive-silicon-dioxide-label`.
    Finally, :ref:`gcmc-silica-label`,
    :ref:`umbrella-sampling-label`,
    are using some more advanced simulation methods
    that are commonly used when studying soft matter systems,
    respectively grand canonical Monte Carlo simulations and 
    a free energy method named umbrella sampling.

Required softwares
==================

..  container:: justify

    The 2 Aug 2023 version of LAMMPS is required to follow the tutorials.
    The other softwares listed here are optional.

LAMMPS (2 Aug 2023)
-------------------

..  container:: justify

    Download and install LAMMPS by following the instructions of the |LAMMPS website installation|.
    Depending on you OS, the procedure may differ.

.. |LAMMPS website installation| raw:: html

   <a href="https://docs.lammps.org/Install.html" target="_blank">LAMMPS website</a>

..  container:: justify

    If you are using *Ubuntu OS*, you can simply execute the
    following command in a terminal:

..  code-block:: bw

   sudo apt-get install lammps

..  container:: justify

    You can verify that LAMMPS is indeed installed on your
    computer by typing in a terminal :

..  code-block:: bw

    lmp

..  container:: justify

    You should see the version of LAMMPS that has been
    installed. On my computer I see

..  code-block:: bw

    LAMMPS (2 Aug 2023)

..  container:: justify

    All the tutorials here were made with the *LAMMPS (2 Aug 2023)*
    version. If you decide to use another LAMMPS version, certain commands
    may not work. In that case, an error message may appear.

VMD (1.9.3)
-----------

..  container:: justify

    In order to visualize the atomic system, the version 1.9.3 of |VMD| will be used.
    If you don't know how to use VMD, you can find some basic instructions here:
    :ref:`vmd-label`. If you prefer, feel free to use an alternative visualization
    software like |Ovito|.
    
.. |VMD| raw:: html

   <a href="https://www.ks.uiuc.edu/Research/vmd" target="_blank">VMD</a>
    
.. |Ovito| raw:: html

   <a href="https://www.ovito.org" target="_blank">Ovito</a>
    
Matplotlib Pyplot
-----------------

..  container:: justify

    In order to plot the results from the simulations,
    you will need a plotting tool. I will use |Matplotlib Pyplot|
    in combination with |lammps_logfile|, a library allowing
    one to read the *log* file produced by LAMMPS.

..  container:: justify

    All the Python scripts used to generate the figure of *LAMMPStutorials*
    can be found in the Github repository.

.. |Matplotlib Pyplot| raw:: html

   <a href="https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html" target="_blank">Matplotlib Pyplot</a>

.. |lammps_logfile| raw:: html

   <a href="https://github.com/henriasv/lammps-logfile" target="_blank">lammps logfile</a>

Text editing software
---------------------

..  container:: justify

    In order to write LAMMPS input file, a text editor is required.
    Any text editor will do, such as |gedit|, |vim|, or |vscode|.
    
.. |gedit| raw:: html

   <a href="https://help.gnome.org/users/gedit/stable/" target="_blank">gedit</a>
    
.. |vim| raw:: html

   <a href="https://www.vim.org/" target="_blank">vim</a>
    
.. |vscode| raw:: html

   <a href="https://code.visualstudio.com/" target="_blank">vscode</a>
    
Find the input scripts
======================

..  container:: justify

    You can access all the input scripts and data files that
    are used in these tutorials from |Github_repository_input| on Github.
    This repository also contains the inputs of every solution to the exercises.

.. |Github_repository_input| raw:: html

    <a href="https://github.com/lammpstutorials/lammpstutorials.github.io/tree/version2.0/docs/inputs" target="_blank">the inputs folder</a>