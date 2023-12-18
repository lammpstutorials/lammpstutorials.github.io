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

Counting the bonds of a CNT
===========================

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

Create a Universe
-----------------

.. container:: justify

    Open a new Jupyter notebook and
    call it *measure_bond_evolution.ipynb*. First,
    let us import both *MDAnalysis*
    and *NumPy* by copying the following
    lines into *measure_bond_evolution.ipynb*.

.. code-block:: python

    import MDAnalysis as mda
    import numpy as np

.. container:: justify

    Then, let us create a MDAnalysis *universe* using
    the LAMMPS data file *cnt_atom.data* as topology,
    and the *lammpstrj* file as trajectory.
    Add the following lines into *measure_bond_evolution.ipynb*:

.. code-block:: python

    path_to_data = "./"
    u = mda.Universe(path_to_data + "cnt_deformed.data",
                     path_to_data + "dump.lammpstrj",
                     topology_format="data", format="lammpsdump",
                     atom_style='id type xs ys zs',
                     guess_bonds=True, vdwradii={'1':1.7})

.. container:: justify

    Since the *.data* file does not contain any bond information
    the original bonds are guessed using the bond guesser
    of MDAnalysis using *guess_bonds=True*.

.. container:: justify

    Note that the bond guesser of MDAnalysis will not update the list of bond
    over time, so we will need to use a few tricks to extract the evolution 
    of the number of bond with time.

.. container:: justify

    Let us create a single atom group
    named *cnt* and containing all the carbon atoms,
    i.e. all the atoms of type 1,
    by adding the following lines into *measure_bond_evolution.ipynb*.

.. code-block:: python

    cnt = u.select_atoms("type 1")

Some basics of MDAnalysis
-------------------------

.. container:: justify

    MDAnalysis allows us to easily access information concerning the simulation, such
    as the number of atoms, or the number of frames in the trajectory:

.. code-block:: python

    print("Number of atoms =", cnt.n_atoms)
    print("Number of frames =", u.trajectory.n_frames)

    Number of atoms = 690
    Number of frames = 286

.. container:: justify

    It is also possible to access the indexes of the atoms that
    are considered as bonded by the bond guesser of MDAnalysis:

.. code-block:: python

    print(cnt.atoms.bonds.indices)

    [[  0   2]
    [  0  23]
    [  0  56]
    (...)
    [686 687]
    [686 689]
    [688 689]]

.. container:: justify

    MDAnalysis also offer the possibility to loop over all the frame of the trajectory using:

.. code-block:: python

    for ts in u.trajectory:
        print(ts.frame)
    
    0
    1
    2
    3
    (...)
    283
    284
    285

.. container:: justify

    The positions of the atoms can also be obtained using:

.. code-block:: python

    u.atoms.positions

    array([[ 75.14728 ,  78.17872 ,  95.61408 ],
    [ 75.33008 ,  77.751114,  93.20232 ],
    [ 75.550476,  77.34152 ,  94.54224 ],
    ...,
    [ 84.66992 ,  82.24888 , 143.84988 ],
    [ 84.66992 ,  82.24888 , 147.60156 ],
    [ 84.85272 ,  81.82128 , 146.26175 ]], dtype=float32)

.. container:: justify

    where the three columns of the array are the *x*,
    *y*,
    and *z* coordinates of the atoms. 

Counting the bonds
------------------

.. container:: justify

    In order to measure the evolution of the number of
    bonds over time, let us loop over the trajectory
    and manually extract the inter-atomic distance over time. 

.. container:: justify

    To to so, for every step of the trajectory, let us
    loop over the indexes of the atoms that were initially
    detected as bonded, and calculate the
    distance between the two atoms, which can be done using:

.. code-block:: python

    for ts in u.trajectory:
        for id1, id2 in cnt.atoms.bonds.indices:
            # detect positions
            pos1 = u.atoms.positions[u.atoms.indices == id1][0]
            pos2 = u.atoms.positions[u.atoms.indices == id2][0]
            r = np.sqrt(np.sum((pos1-pos2)**2))

.. container:: justify

    Then, let us assume that if :math:`r` is larger that a 
    certain cut-off value of, let's say, :math:`1.8\,Å`,
    the bond is broken:

.. code-block:: python

    for ts in u.trajectory:
        for id1, id2 in cnt.atoms.bonds.indices:
            pos1 = u.atoms.positions[u.atoms.indices == id1][0]
            pos2 = u.atoms.positions[u.atoms.indices == id2][0]
            r = np.sqrt(np.sum((pos1-pos2)**2))
            if r < 1.8:
                print("the bond has a length", r, "Å")
            else:
                print("the bond is broken")
    
.. container:: justify

    Finally, let us store both mean length of bonds
    and total number of bond in lists.  

.. code-block:: python

    lbond_vs_frame = []
    nbond_vs_frame = []
    for ts in u.trajectory:
        frame = ts.frame
        all_bonds_ts = [] # temporary list to store bond length
        for id1, id2 in cnt.atoms.bonds.indices:
            pos1 = u.atoms.positions[u.atoms.indices == id1]
            pos2 = u.atoms.positions[u.atoms.indices == id2]
            r = np.sqrt(np.sum((pos1-pos2)**2))
            if r < 1.8:
                all_bonds_ts.append(r)
        mean_length_bonds = np.mean(all_bonds_ts)
        number_of_bond = len(all_bonds_ts)/2 # divide by 2 to avoid counting twice
        lbond_vs_frame.append([frame, mean_length_bonds]) 
        nbond_vs_frame.append([frame, number_of_bond])

.. container:: justify

    The data can then be saved to files:

.. code-block:: python

    np.savetxt("number_bond_vs_time.dat", nbond_vs_time)
    np.savetxt("length_bond_vs_time.dat", lbond_vs_time)

.. figure:: ../figures/mdanalysis/bond-dark.png
    :alt: plot of the bond length and distance versus time
    :class: only-dark

.. figure:: ../figures/mdanalysis/bond-light.png
    :alt: plot of the bond length and distance versus time
    :class: only-light

.. container:: figurelegend

    Figure: Evolution of the average bond length (a) and bond number (b) as a function of time.

Bond length distributions
-------------------------

.. container:: justify

    Using a similar script,
    let us extract the bond length distribution
    at the beginning of the simulation (let us say the 20 first frame),
    as well as near the maximum deformation of the CNT:

.. code-block:: python

    for ts in u.trajectory:
        all_bonds_ts = []
        for id1, id2 in cnt.atoms.bonds.indices:
            pos1 = u.atoms.positions[u.atoms.indices == id1]
            pos2 = u.atoms.positions[u.atoms.indices == id2]
            r = np.sqrt(np.sum((pos1-pos2)**2))
            if r < 1.8:
                all_bonds_ts.append(r)
        if frame > 0: # ignore the first frame
            histo, r_val = np.histogram(all_bonds_ts, bins=50, range=(1.3, 1.65))
            r_val = (r_val[1:]+r_val[:-1])/2
            bond_length_distributions.append(np.vstack([r_val, histo]))

.. figure:: ../figures/mdanalysis/bond-distribution-dark.png
    :alt: plot of the bond distribution
    :class: only-dark

.. figure:: ../figures/mdanalysis/bond-distribution-light.png
    :alt: plot of the bond distribution
    :class: only-light

.. container:: figurelegend

    Figure: Distribution in bond length near the start of the simulation,
    as well as near the maximum deformation of the CNT.