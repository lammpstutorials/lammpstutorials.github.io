[![DOI](https://zenodo.org/badge/737122279.svg)](https://zenodo.org/doi/10.5281/zenodo.13341527)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Matplotlib Pyplot functions

My own Python class for making [Pyplot](https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html)
figures compatible with Latex documents. These functions are used in my most recent publications, as well as in:
- [LAMMPS tutorials](https://lammpstutorials.github.io)
- [NMRforMD](https://nmrformd.readthedocs.io)

## Examples

### Simple bi-panel figure

![illustration](examples/bi-panel.png#gh-light-mode-only)

![illustration](examples/bi-panel-dm.png#gh-dark-mode-only)

The Python code to generate the figure:

```python
    # Initialise figure
    myplt = PltTools()
    myplt.prepare_figure(fig_size = (18,12), dark_mode = dark_mode,
                         transparency = True, n_colone = 2)
    # Panel 1
    myplt.add_panel(panel_position = [0, 0])
    myplt.add_plot(x = x1, y = y1/2, marker = "o", data_color = 2,
                   markersize = 6, data_label = r"$\textrm{data}~\alpha$")
    myplt.add_plot(x = x3, y = y3, marker = "o", data_color = 0, markersize = 10,
                   data_label = r"$\textrm{data}~\beta$")
    myplt.add_plot(x = x2, y = y2, marker = "-", data_color = "autogray",
                   data_label = r"$\textrm{linear fit}$")
    myplt.complete_panel(xlabel = r"$x ~ \textrm{[\AA{}]}$", ylabel = r"$t ~ \textrm{[ps]}$",
                         panel_title = r"$\textrm{Time vs distance}$", xpad = 10)
    myplt.set_boundaries(y_boundaries = (0, 5))
    # Panel 2
    myplt.add_panel(type = 'loglog', panel_position = [0, 1])
    myplt.add_plot(x = x4, y = y4, marker = "o", data_color = 0,
                   markersize = 16)
    myplt.add_plot(x = x5, y = y5, marker = "p", data_color = 0,
                   markersize = 16)
    myplt.add_plot(x = x6, y = y6, marker = "^", data_color = 0,
                   markersize = 16)
    myplt.complete_panel(xlabel = r"$E ~ \textrm{[J]}$", ylabel = r"$N_\textrm{p} ~ [0]$",
                         panel_title = r"$\textrm{Population size vs energy}$", xpad = 10)
    myplt.set_boundaries(x_boundaries = (0.01, 1000), y_boundaries = (0.01, 1000))
    # Finish figure
    myplt.add_subplotlabels(type_label_panel = "a")
    myplt.save_figure(filename = "bi-panel", saving_path = 'examples/', show = False)
```

See also the [examples](examples/examplesmv.ipynb) notebook.

## Light mode vs dark mode

- For web integration, use either dark or light mode with a transparent background, as done for instance in [lammps tutorials](https://lammpstutorials.github.io).
- For scientific publication, use the light mode without transparent background.
