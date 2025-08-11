.. _using-lammps-gui-label:

Using LAMMPS--GUI
*****************

.. admonition:: Note
    :class: non-title-info

    For simplicity, these tutorials reference keyboard shortcuts
    based on the assignments for Linux and Windows. macOS users should
    use the ``Command`` key in place of the
    ``Ctrl`` key when using keyboard shortcuts.

Installation
============

Precompiled versions of LAMMPS--GUI are available for Linux, macOS,
and Windows on the LAMMPS GitHub Release
page.  The Linux version is provided in two
formats: as compressed tar archive (.tar.gz) and as a |flatpak|
bundle.  The macOS version is distributed as a
.dmg installer image, while the Windows version comes as an executable
installer package.

.. |flatpak| raw:: html

   <a href="https://flatpak.org" target="_blank">Flatpak</a>

Installing the Linux .tar.gz Package
------------------------------------

Download the archive (e.g., LAMMPS-Linux-x86_64-GUI-22Jul2025.tar.gz)
and unpack it.  This will create a folder named LAMMPS--GUI containing the
included commands, which can be launched directly using ``./lammps-gui`` or
``./lmp``, for example.  Adding this folder to the PATH environment
variable will make these commands accessible from everywhere, without the
need for the ``./`` prefix.

Installing the Linux Flatpak Bundle
-----------------------------------

You have to have Flatpak support installed on Linux machine to be able
to use the Flatpak bundle.  Download the bundle file
(e.g., LAMMPS-Linux-x86_64-GUI-22Jul2025.flatpak) and then
install it using the following command:

.. code-block:: bash

    flatpak install --user LAMMPS-Linux-x86_64-GUI-22Jul2025.flatpak

This will integrate LAMMPS--GUI into your desktop environment
(e.g., GNOME, KDE, XFCE) where it should appear in the ``Applications``
menu under ``Science``.  Additionally, the ``.lmp`` file extension will be
registered to launch LAMMPS--GUI when opening a file with this
extension in the desktop's file manager.

You can also launch LAMMPS--GUI from the command-line using the following command:

.. code-block:: bash

    flatpak run org.lammps.lammps-gui

Similarly, for launching the LAMMPS command-line executable, use:

.. code-block:: bash
    
    flatpak run --command=lmp org.lammps.lammps-gui -in in.lmp

Installing the macOS Application Bundle
---------------------------------------

After downloading the macOS app bundle image file
(e.g., LAMMPS-macOS-multiarch-GUI-22Jul2025.dmg), double-click
on it.  In the dialog that opens drag the LAMMPS--GUI app bundle into
the Applications folder.  To enable command-line access, follow the
instructions in the **README.txt** file.  These macOS app-bundles contain
native executables for both, Intel and Apple CPUs.

After installation, you can launch LAMMPS--GUI from the Applications
folder.  Additionally, you can drag an input file onto the app or open
files with the ``.lmp`` extension.  Note that the LAMMPS--GUI app bundle is
currently not cryptographically signed, so macOS may initially prevent
it from launching.  If this happens, you need to adjust the settings in
the ``Security & Privacy`` system preferences dialog to allow access.

Installing the Windows package
------------------------------

Download the LAMMPS--GUI installer for Windows
(e.g., LAMMPS-Win10-64bit-GUI-22Jul2025.exe).  Windows may warn
you that the file is from an unknown developer and was downloaded from
the internet.  This happens because neither the installer nor the
LAMMPS--GUI application (or any other included applications) have been
cryptographically signed.  You will need to choose to keep the file, and
when launching the installer, confirm that you want to run it despite
the warning.

After installation, a new entry should appear in the Start menu.
Additionally, the ``.lmp`` file extension should be registered with
Windows File Explorer to open LAMMPS--GUI when opening a file with the
``.lmp`` extension.  The ``lammps-gui`` and ``lmp`` commands should also
be available in the command-line.

Opening, Editing, and Saving Files
----------------------------------

LAMMPS--GUI can be launched from the command-line, as explained above, where you
can either launch it without arguments or provide one file name as an argument.  All
other arguments will be ignored.  For example:

.. code-block:: bash

    lammps-gui input.lmp

Files can also be opened from the ``File`` menu.  You can select a
file through a dialog and then open it.  Additionally, a history of
the last five opened files is maintained, with entries to open them directly.
Finally, the ``Ctrl-O`` keyboard shortcut can also be used to open a file.

When integrated into a desktop environment, it is also possible to open
files with a ``.lmp`` extension or use drag-and-drop.

For the most part, the editor window behaves like other graphical
editors.  You can enter, delete, or copy and paste text.   When entering
text, a pop-up window will appear with possible completions after typing
the first two characters of the first word in a line.  You can
navigate the highlighted options using the up and down arrow keys, and select a
completion by pressing the Enter key or using the mouse.  You can also continue
typing, and the selection in the pop-up will be refined.  For some
commands, there will be completion pop-ups for their
keywords or when a filename is expected, in which case,
the pop-up will list files in the current folder.

As soon as LAMMPS--GUI recognizes a command, it applies syntax
highlighting according to built-in categories.  This can help
detect typos, since those may cause LAMMPS--GUI not to
recognize the syntax and thus not apply or partially apply
the syntax highlighting.  When you press the ``Tab`` key, the line will be
reformatted.  Consistent formatting can improve the readability of
input files, especially long and complex ones.

If the file in the editor has unsaved changes, the word
*modified* will appear in the window title.  The current input
buffer can be saved by selecting ``Save`` or ``Save As...`` from the
``File`` menu.  You can also click the ``Save`` icon on the left side
of the status bar, or use the ``Ctrl-S`` keyboard shortcut.

.. admonition:: Note
    :class: non-title-info

    When LAMMMPS--GUI opens a file, it will *switch* the working directory
    to the folder that contains the input file.  The same happens when saving to
    a different folder than the current working directory.  The current working
    directory can be seen in the status bar at the bottom right.  This is important
    to note because LAMMPS input files often require additional files for reading and may
    write output files (such as images, trajectory dumps, or averaged data files),
    which are typically expected to be in the same folder as the input file.
