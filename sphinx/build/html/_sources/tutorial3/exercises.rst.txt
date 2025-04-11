Going further with exercises
============================

Extract the radial distribution function
----------------------------------------

Extract the radial distribution functions (RDF or :math:`g(r)`)
between the oxygen atom of the water molecules
and the oxygen atom from the PEG molecule. Compare the rdf
before and after the force is applied to the PEG.

.. figure:: figures/RDF-dark.png
    :alt: RDF g(r) for water and peg
    :class: only-dark

.. figure:: figures/RDF-light.png
    :alt: RDF g(r) for water and peg
    :class: only-light

.. container:: figurelegend

    Figure: Radial distribution function between the oxygen atoms 
    of water, as well as between the oxygen atoms of water and the 
    oxygen atoms of the PEG molecule.  

Note the difference in the structure of the water before and after
the PEG molecule is stretched. This effect is described in
the 2017 publication by Liese et al. :cite:`liese2017hydration`.

Add salt to the system
----------------------

Realistic systems usually contain ions. Let us add some :math:`\text{Na}^+` and 
:math:`\text{Cl}^-` ions to our current PEG-water system.

Add some :math:`\text{Na}^+` and 
:math:`\text{Cl}^-` ions to the mixture using the method
of your choice. :math:`\text{Na}^+` ions are 
characterised by their mass :math:`m = 22.98\,\text{g/mol}`,
their charge :math:`q = +1\,e`, and Lennard-Jones
parameters, :math:`\epsilon = 0.0469\,\text{kcal/mol}`
and :math:`\sigma = 0.243\,\text{nm}`,
and :math:`\text{Cl}^-` ions by their
mass :math:`m = 35.453\,\text{g/mol}`,
charge :math:`q = -1\,e` and Lennard-Jones
parameters, :math:`\epsilon = 0.15\,\text{kcal/mol}`,
and :math:`\sigma = 0.4045\,\text{nm}`.

.. figure:: figures/salt-exercise-dark.png
    :alt: PEG in a NaCl solution
    :class: only-dark

.. figure:: figures/salt-exercise-light.png
    :alt: PEG in a NaCl solution
    :class: only-light

.. container:: figurelegend

    Figure: A PEG molecule in the electrolyte with :math:`\text{Na}^+` ions in 
    purple and :math:`\text{Cl}^-` ions in cyan.

Evaluate the deformation of the PEG
-----------------------------------

Once the PEG is fully stretched, its structure differs from the
unstretched case. The deformation can be probed by extracting the typical
intra-molecular parameters, such as the typical angles of the dihedrals.

Extract the histograms of the angular distribution of the PEG dihedrals
in the absence and the presence of stretching.

.. figure:: figures/dihedral_angle-dark.png
    :alt: PEG in a NaCl solution
    :class: only-dark

.. figure:: figures/dihedral_angle-light.png
    :alt: PEG in a NaCl solution
    :class: only-light

.. container:: figurelegend

    Figure: Probability distribution for the dihedral angle :math:`\phi`, for a stretched
    and for an unstretched PEG molecule.
