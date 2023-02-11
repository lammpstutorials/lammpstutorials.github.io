.. contribute-label:

Wanna contribute?
*****************

Your feedback is always useful. A particular thanks to Andrea Corradini who spontaneously send me a very detailed report of the 
:ref:`umbrella-sampling-label` tutorial, and signaled some inconsistencies.

You can also suggest idea for new a tutorial, or a new section to an existing tutorial, just keep in mind that:

- the content must be of interest for everyone, not just for the specific issue you are trying to solve,
- the tutorials must remain relatively simple. 

As a tester
===========

Report broken link and typo by `email`_ or by posting a new issue on `github`_.

.. _email: simon.gravelle@live.fr
.. _github: https://github.com/lammpstutorials/lammpstutorials.github.io/issues

As a writer
===========

Propose a new tutorial or a modification to an existing tutorial.
To do so, fork the lammpstutorials repository on github, make your changes,
and submit a `pull request`_.

.. _pull request: https://github.com/lammpstutorials/lammpstutorials.github.io/pulls

Build lammpstutorial locally
============================

*lammpstutorial* can be build locally on your computer using sphinx. Fork
the repository by typing in a terminal:

..  code-block:: bash

    git clone https://github.com/lammpstutorials/lammpstutorials.github.io.git

Then, go to docs/sphinx/, and execute:

..  code-block:: bash

    bash build.sh

The docs/index.html can be opened with a web browser.
The tutorials are located in docs/sphinx/source/tutorials and written in rst format. 
For each tutorials, the corresponding lammps input files are located here: docs/inputs.