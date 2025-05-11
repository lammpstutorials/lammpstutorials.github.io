.. _running-lammps-label:

Running LAMMPS
==============

From within the LAMMPS--GUI main window, LAMMPS can be started either from
the ``Run`` menu by selecting the ``Run LAMMPS from Editor Buffer`` entry,
using the keyboard shortcut Ctrl-Enter (Command-Enter on macOS), or by clicking the
green ``Run`` button in the status bar.  While LAMMPS is running, a message on
the left side indicates that LAMMPS is active, along with the number of active threads.
On the right side, a progress bar is displayed, showing the estimated progress
of the current ``run`` or ``minimize`` command.

Creating Snapshot Images
------------------------

Open the ``Image Viewer`` using either the ``Create Image`` option
from the ``Run`` menu, the ``Ctrl-I`` keyboard shortcut,
or click on the (right) palette button in the status bar.  The image
can be saved using the ``Save As...`` option from the ``File`` menu.

The Output Window
-----------------

By default, when starting a run, the ``Output`` window opens to display the screen
output of the running LAMMPS calculation.  The text in the Output window is
read-only and cannot be modified, but keyboard shortcuts for selecting and
copying all or part of the text can be used to transfer it to another program:
The keyboard shortcut ``Ctrl-S`` (or ``Command-S`` on ``macOS``) can
be used to save the Output buffer to a file.  Additionally, the ``Select All``
and ``Copy`` functions, along with a ``Save Log to File`` option, are available
through the context menu, which can be accessed by right-clicking within the text area of the
``Output`` window.

The Charts Window
-----------------

By default, when starting a run, a ``Charts`` window opens to display
a plot of the thermodynamic output from the LAMMPS calculation.  From the ``File``
menu in the top-left corner, you can save an image of the
currently displayed plot or export the data in various formats:
plain text columns (for use with plotting tools like Gnuplot or XmGrace),
CSV data (suitable for processing in Microsoft Excel, LibreOffice Calc,
or Python with Pandas), or YAML (which can be imported into Python using PyYAML or Pandas).
You can use the mouse to zoom in on the graph by holding the left button and dragging
to select an area.  To zoom out, right-click anywhere on the graph.  You can reset the view
by clicking the ``lens`` button located next to the data drop-down menu.

Preferences
-----------

The Preferences dialog allows customization of the behavior and appearance of
LAMMPS--GUI.  Among other options:

- In the ``General Settings`` tab, the ``Data update interval`` setting
  allows you to define the time interval, in milliseconds, between data updates during
  a LAMMPS run.  By default, the data for the ``Charts`` and ``Output``
  windows is updated every 10 milliseconds.  Set this to 100 milliseconds or more
  if LAMMPS--GUI consumes too many resources during a run.  The ``Charts update interval``
  controls the time interval between redrawing the plots in the ``Charts`` window, in milliseconds.
- The ``Accelerators`` tab enables you to select an accelerator package
  for LAMMPS to use.  Only settings supported by the LAMMPS library and local hardware
  are available.  The ``Number of threads`` field allows you to set the maximum
  number of threads for accelerator packages that utilize threading.
- The ``Editor Settings`` tab allows you to adjust the settings of the editor
  window.  Select the ``Auto-save on Run and Quit`` option to automatically save changes
  made to the ``.lmp`` file upon closing LAMMPS--GUI.

See |LAMMPSGUIDOC| for a full list of options.

.. |LAMMPSGUIDOC| raw:: html

   <a href="https://docs.lammps.org/Howto_lammps_gui.html" target="_blank">How-to LAMMPS--GUI</a>