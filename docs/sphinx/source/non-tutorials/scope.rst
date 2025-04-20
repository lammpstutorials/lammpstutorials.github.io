.. _scope-label:

Scope
*****

This set of tutorials consists of eight tutorials arranged in order of
increasing difficulty.  The novelties associated with each tutorial are
briefly described below.

In :ref:`lennard-jones-label`, the structure of LAMMPS
input files is illustrated through the creation of a simple atomic
Lennard-Jones fluid system.  Basic LAMMPS commands are used to set up
interactions between atoms, perform an energy minimization, and finally
run a simple MD simulation in the microcanonical (NVE) and canonical (NVT)
ensembles.

In :ref:`carbon-nanotube-label`, a more complex system
is introduced in which atoms are connected by bonds: a small carbon
nanotube.  The use of both classical and reactive force fields (here,
OPLS-AA :cite:`jorgensenDevelopmentTestingOPLS1996` and
AIREBO :cite:`stuart2000reactive`, respectively) is illustrated.  An
external deformation is applied to the CNT, and its deformation is
measured.  This tutorial also demonstrates the use of an external tool
to visualize breaking bonds, and show the possibility to import
LAMMPS-generated YAML log files into Python.

In :ref:`all-atoms-label`, two components\textemdash liquid water
(flexible three-point model) and a polymer molecule\textemdash are merged and
equilibrated.  A long-range solver is used to handle the electrostatic
interactions accurately, and the system is equilibrated in the
isothermal-isobaric (NPT) ensemble; then, a stretching force is applied
to the polymer.  Through this relatively complex solvated polymer
system, the tutorial demonstrates how to use type labels to make
molecule files more generic and easier to manage :cite:`gissinger2024type`.

In :ref:`sheared-confined-label`, an electrolyte is
confined between two walls, illustrating the specifics of simulating
systems with fluid-solid interfaces.  With the rigid four-point
TIP4P/2005 :cite:`abascal2005general` water model, this tutorial uses a
more complex water model than :ref:`all-atoms-label`.  A
non-equilibrium MD is performed by imposing shear on the fluid through
moving the walls, and the fluid velocity profile is extracted.

In :ref:`reactive-silicon-dioxide-label`, the ReaxFF
reactive force field is used, specifically designed to simulate chemical
reactions by dynamically adjusting atomic interactions
:cite:`van2001reaxff`.  ReaxFF includes charge equilibration (QEq), a
method that allows the partial charges of atoms to adjust according to
their local environment.

In :ref:`gcmc-silica-label`, a Monte Carlo simulation in
the grand canonical ensemble is implemented to demonstrate how LAMMPS
can be used to simulate an open system that exchanges particles with a
reservoir.

In :ref:`umbrella-sampling-label`, an advanced free
energy method called umbrella sampling is implemented.  By calculating
an energy barrier, this tutorial describes a protocol
for addressing energy landscapes that are difficult to sample using
classical MD or MC methods.

In :ref:`bond-react-label`, a CNT embedded in
nylon-6,6 polymer melt is simulated.  The
REACTER protocol is used to model the polymerization of nylon, and the formation
of water molecules is tracked over time~\cite{gissinger2020reacter}.