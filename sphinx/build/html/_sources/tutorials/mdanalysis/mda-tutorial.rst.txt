.. _mda-label:

MDAnalysis tutorial
*******************

.. container:: hatnote

    Perform post-mortem analysis using Python

.. figure:: ../figures/level1/breaking-a-carbon-nanotube/CNT_dark.webp
    :alt: carbon nanotube image in vacuum
    :height: 250
    :align: right
    :class: only-dark

.. figure:: ../figures/level1/breaking-a-carbon-nanotube/CNT_light.webp
    :alt: carbon nanotube image in vacuum
    :height: 250
    :align: right
    :class: only-light

.. container:: justify

    There are two main ways to analyze data from a MD simulation:
    (1) on-the-fly analysis, like is done for instance using *fix ave/time*,
    and (2) post-mortem analysis. Post-mortem analysis can be performed using
    the atom coordinates saved in the *lammpstrj* file. 

.. container:: justify

    The goal of this extra tutorial is to provide some tips
    to import *lammpstrj*  trajectory file into Python.

.. include:: ../../non-tutorials/needhelp.rst

Counting bonds
==============

.. container:: justify

    Here, we re-use the trajectory generated
    during the second part *Breakable bonds*
    of the :ref:`carbon-nanotube-label` tutorial.
    It is recommended that you follow this tutorial
    first, but you can also directly download the |dump_cnt|
    file and the |data_cnt| file and continue with this MDA tutorial.

.. |dump_cnt| raw:: html

   <a href="../../../../../inputs/level1/breaking-a-carbon-nanotube/breakable-bonds/dump.lammpstrj" target="_blank">dump</a>

.. |data_cnt| raw:: html

   <a href="../../../../../inputs/level1/breaking-a-carbon-nanotube/breakable-bonds/cnt_atom.data" target="_blank">data</a>

.. container:: justify

    Open a new Jupyter notebook and
    call it *measure_bond_evolution.ipynb*. First, let us import both *MDAnalysis*
    and *NumPy* by copying the following lines into *measure_bond_evolution.ipynb*.

.. code-block:: python

    import MDAnalysis as mda
    import numpy as np

.. container:: justify

    Then, let us create a *MDAnalysis* *universe* using the LAMMPS
    data file *cnt_atom.data* as topology,
    and the *lammpstrj* file as trajectory. Let us guess the
    original bonds using the bond guesser of MDAnalysis (*guess_bonds=True*).
    Let us also create a single atom group named *cnt* and containing all the carbon atoms.
    Add the following lines into *measure_bond_evolution.ipynb*:

 .. code-block:: python

    path_to_data = "./"
    u = mda.Universe(path_to_data + "cnt_deformed.data",
                     path_to_data + "dump.lammpstrj",
                     topology_format="data", format="lammpsdump",
                     atom_style='id type xs ys zs',
                     guess_bonds=True, vdwradii={'1':1.7})
    cnt = u.select_atoms("type 1")

.. container:: justify

    Note that the bond guesser of MDAnalysis will not update the list of bond
    over time, so we will need to use a few tricks to extract the evolution 
    of the number of bond with time.

.. container:: justify

    Then, let us loop over the trajectory and extract the bonds average length
    and total number over time. Add the
    following lines into *measure_bond_evolution.ipynb*:

.. code-block:: python

        nbond_vs_time = []
        lbond_vs_time = []
        # loop over trajectory
        for ts in u.trajectory:
            # sabe the bond of the timestep ts in a list
            all_bonds_ts = []
            # loop over all initially detected bond
            for id1, id2 in cnt.atoms.bonds.indices:
                # detect positions
                pos1 = u.atoms.positions[u.atoms.indices == id1]
                pos2 = u.atoms.positions[u.atoms.indices == id2]
                d = pos1-pos2
                r = np.sqrt(d[:, 0]**2 + d[:, 1]**2 + d[:, 2]**2)
                if r < 1.8: # assume that bond longer than 1.8 angstroms are broken
                    all_bonds_ts.append(r)
            lbond_vs_time.append([ts.time*5000*0.0005, np.mean(all_bonds_ts)]) 
            nbond_vs_time.append([ts.time*5000*0.0005, len(all_bonds_ts)/2]) # divide by 2 to avoid counting twice
        nbond_vs_time = np.array(nbond_vs_time)
        lbond_vs_time = np.array(lbond_vs_time)
        # finally let us save the data to tex files:
        np.savetxt("number_bond_vs_time.dat", nbond_vs_time)
        np.savetxt("length_bond_vs_time.dat", lbond_vs_time)

.. container:: justify

    The array *nbond_vs_time* contains the number of bond as a function of time, and 
    *lbond_vs_time* the bond length. Let us plot both of them:

.. figure:: ../figures/mdanalysis/bond-dark.png
    :alt: plot of the bond length and distance versus time
    :class: only-dark

.. figure:: ../figures/mdanalysis/bond-light.png
    :alt: plot of the bond length and distance versus time
    :class: only-light

.. container:: figurelegend

    Figure: Evolution of the bond length (a) and bond number (b) as a function of time.