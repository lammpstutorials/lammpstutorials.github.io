.. _prerequisites-label:

Prerequisites
*************

Background knowledge
====================

This set of tutorials assumes no prior knowledge of the LAMMPS software
itself.  To complete the tutorials, a text editor and a suitable LAMMPS
executable are required.  We use |lammpsguidocs|
here, as it offers features that make it particularly convenient for
tutorials, but other console or graphical text editors, such as GNU
nano, vi/vim, Emacs, Notepad, Gedit, and Visual Studio Code can also be
used.  LAMMPS can be executed either directly from LAMMPS--GUI
(:ref:`using-lammps-gui-label`) or from a command prompt
(:ref:`command-line-label`), the latter of which requires some familiarity
with executing commands from a terminal or command-line prompt.

In addition, prior knowledge of the theoretical basics of molecular
simulations and statistical physics is highly beneficial.  Users may
refer to textbooks such as *Understanding Molecular Simulation* by
Daan Frenkel and Berend Smit :cite:`frenkel2023understanding`, as well as
*Computer Simulation of Liquids* by Michael Allen and Dominic
Tildesley :cite:`allen2017computer`.  To better understand
the fundamental concepts behind the soft matter systems simulated in these
tutorials, users can also refer to *Basic Concepts for Simple and Complex Liquids*
by Jean-Louis Barrat and Jean-Pierre Hansen
:cite:`barrat2003basic`, as well as
*Theory of Simple Liquids: with Applications to Soft Matter*
by Jean-Pierre Hansen and Ian Ranald McDonald :cite:`hansen2013theory`. 
For more resources, the |sklogwiki_main_page| platform provies a wide range of information
on statistical mechanics and molecular simulations.

.. |sklogwiki_main_page| raw:: html

    <a href="http://www.sklogwiki.org/SklogWiki/index.php/Main_Page" target="_blank">SklogWiki</a>

Software/system requirements
============================

The LAMMPS stable release version 22Jul2025
and the matching LAMMPS--GUI software version 1.6.12 are required to
follow the tutorials, as they include features that were first
introduced in these versions.  For Linux (x86_64 CPU), macOS (BigSur or
later), and Windows (10 and 11) you can download a precompiled LAMMPS
package from the LAMMPS |LAMMPSrelease| page on
GitHub.  Select a package with ``GUI`` in the
file name, which includes both, LAMMPS--GUI and the LAMMPS command-line
executable.  These precompiled packages are designed to be portable, and
therefore omit support for parallel execution with MPI.  Instructions
for installing LAMMPS--GUI and using its most relevant features for the
tutorials are provided in :ref:`using-lammps-gui-label`.

.. |LAMMPSrelease| raw:: html

    <a href="https://github.com/lammps/lammps/releases" target="_blank">release</a>

LAMMPS versions are generally backward compatible, meaning that older
input files typically work the same with newer versions of LAMMPS.
However, forward compatibility is not as strong, so input files written
for a newer version may not always work with older versions.  As a
result, it is usually possible to follow this tutorial with more recent
releases of LAMMPS--GUI and LAMMPS; older versions may require some
(minor) adjustments.  These tutorials will be periodically updated to
ensure compatibility and benefit from new features in the latest stable
version of LAMMPS.

For some tutorials, external tools are required for plotting and
visualization, as the corresponding functionality in LAMMPS--GUI is
limited.  Suitable tools for plotting include Python with
Pandas/Matplotlib :cite:`van1995python` :cite:`hunter2007Matplotlib`, XmGrace,
Gnuplot, Microsoft Excel, or LibreOffice Calc.  For visualization,
suitable tools include |vmd| :cite:`humphrey1996vmd` and
|ovito| :cite:`stukowski2009visualization`.

.. |ovito| raw:: html

    <a href="https://ovito.org" target="_blank">OVITO</a>

.. |vmd| raw:: html

    <a href="https://www.ks.uiuc.edu/Research/vmd" target="_blank">VMD</a>

About LAMMPS--GUI
=================

LAMMPS--GUI is a graphical text editor, enhanced for editing LAMMPS
input files and linked to the LAMMPS library, allowing it to run LAMMPS
directly.  The text editor functions similarly to other graphical
editors, such as Notepad or Gedit, but offers the following enhancements
specifically for LAMMPS:

- Wizard dialogs to set up these tutorials
- Auto-completion of LAMMPS commands and options
- Context-sensitive online help
- Syntax highlighting for LAMMPS input files
- Syntax-aware line indentation
- Visualization using LAMMPS' built-in renderer
- Start and stop simulations via mouse or keyboard
- Monitoring of simulation progress
- Dynamic capture of LAMMPS output in a text window
- Automatic plotting of thermodynamic data during runs
- Capture of ``dump image`` outputs for animations
- Export of thermodynamic data for external plotting
- Inspection of binary restart files

:ref:`using-lammps-gui-label` contains basic instructions for installation and using LAMMPS--GUI with
the tutorials presented here.  A complete description of all LAMMPS--GUI
features can be found in the LAMMPS manual (see |lammpsguidocs|).

.. |lammpsguidocs| raw:: html

    <a href="https://docs.lammps.org/stable/Howto_lammps_gui.html" target="_blank">LAMMPS--GUI</a>

Content and citation
====================

All files and inputs required to follow the tutorials are available from a
dedicated GitHub organization account, |lammpstutorials_organization|.
If you find these tutorials useful, you can
cite *A Set of Tutorials for the LAMMPS Simulation Package* by Simon Gravelle,
Jacob R. Gissinger, and Axel Kohlmeyer (2025) :cite:`gravelle2025tutorials`. You
can access the full paper on |gravelle2025tutorials_arXiv|.

.. |lammpstutorials_organization| raw:: html

    <a href="https://github.com/lammpstutorials" target="_blank">LAMMPStutorials</a>

.. |gravelle2025tutorials_arXiv| raw:: html

    <a href="https://doi.org/10.48550/arXiv.2503.14020" target="_blank">arXiv</a>
