# Matplotlib Pyplot functions

Personal functions for making [Pyplot](https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html) figures compatible with Latex documents.
These functions are used for several of my personal projects, including [LAMMPS tutorials](https://lammpstutorials.github.io) and 
[NMRforMD](https://nmrformd.readthedocs.io).

## Light mode vs dark mode

Most of my [projects](https://lammpstutorials.github.io) have a dark and a light modes, 
therefore I systematically generate figures for both modes, with 
transparent background. For a scientific publication, I use the light mode
without transparent background.

## Example

See the [examples](examples.ipynb) notebook for the commands behind these figures.

### Single panel with legend

![illustration](examples/example-1-light.png#gh-light-mode-only)

![illustration](examples/example-1-dark.png#gh-dark-mode-only)

### Triple panels with two different widths and labels

![illustration](examples/example-2-light.png#gh-light-mode-only)

![illustration](examples/example-2-dark.png#gh-dark-mode-only)
