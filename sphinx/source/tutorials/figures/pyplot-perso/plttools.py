import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.ticker import AutoMinorLocator

sys.path.append("colors/")
from colorseries import colorserie1, gray_for_dm, gray_for_lm

class PltTools():
    """
    fig_size must be chosen as a tuple, e.g. (18,6)
    If dark_mode, then a dark background is set.
    """
    def __init__(self,
                 fig_size = None,
                 dark_mode = False,
                 transparency = False,
                 use_serif = True,
                 tex_font = True,
                 panel_position = None,
                 x = None,
                 y = None,
                 marker = None,
                 type = "plot",
                 markersize = 12,
                 data_linewidth = 4,
                 data_color = None,
                 axis_color = None,
                 data_label = None,
                 panel_label = None,
                 shift_label_panel = 0.2,
                 sshift_label_panel = None,
                 type_label_panel = "a",
                 open_symbols = False,
                 markeredgewidth = 0,
                 n_colone = 1,
                 n_line = 1,
                 xlabel = None,
                 ylabel = None,
                 cancel_x = False,
                 cancel_y = False,
                 fontsize = 'default',
                 tickwidth1 = 2.5,
                 tickwidth2 = 2,
                 legend = True,
                 ncol_legend = 1,
                 locator_x = 'auto',
                 locator_y = 'auto',
                 panel_title = None,
                 panel_linewidth = 2.5,
                 xpad = 10,
                 ypad = 6,
                 x_boundaries=None,
                 x_ticks=None,
                 y_boundaries=None,
                 y_ticks=None,
                 saving_path = None,
                 filename = None,
                 *args,
                 **kwargs
                 ):
        super().__init__(*args, **kwargs)
        self.fig_size = fig_size
        self.dark_mode = dark_mode
        self.transparency = transparency
        self.use_serif = use_serif
        self.tex_font = tex_font
        self.x = x
        self.y = y
        self.panel_position = panel_position
        self.type = type
        self.marker = marker
        self.markersize = markersize
        self.data_linewidth = data_linewidth
        self.data_color = data_color
        self.axis_color = axis_color
        self.data_label = data_label
        self.panel_label = panel_label
        self.open_symbols = open_symbols
        self.markeredgewidth = markeredgewidth
        self.shift_label_panel = shift_label_panel
        self.sshift_label_panel = sshift_label_panel
        self.type_label_panel = type_label_panel
        self.n_colone = n_colone
        self.n_line = n_line
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cancel_x = cancel_x
        self.cancel_y = cancel_y
        self.fontsize = fontsize
        self.tickwidth1 = tickwidth1
        self.tickwidth2 = tickwidth2
        self.legend = legend
        self.ncol_legend = ncol_legend
        self.locator_x = locator_x
        self.locator_y = locator_y
        self.panel_title = panel_title
        self.panel_linewidth = panel_linewidth
        self.xpad = xpad
        self.ypad = ypad
        self.x_boundaries=x_boundaries
        self.x_ticks=x_ticks
        self.y_boundaries=y_boundaries
        self.y_ticks=y_ticks
        self.saving_path = saving_path
        self.filename = filename

        if self.fontsize == 'default':
            self.fontsize = 34
            self.font = {'family': 'sans',
                         'color':  'black',
                         'weight': 'normal',
                         'size': self.fontsize}
        
    def update_parameters(self, **args):    
        for arg, value in args.items():
            setattr(self, arg, value)
        
        if self.axis_color is None:
            if self.dark_mode:
                self.axis_color = gray_for_dm
            else:
                self.axis_color = gray_for_lm

    def reset_parameters(self):
        """Some parameters value need to be forgotten."""
        self.open_symbols = False
        self.data_label = False
        self.marker = None
        self.type = "plot"
        self.markersize = 12
        self.data_linewidth = 4
        self.data_color = None

    def prepare_figure(self, **args):

        self.update_parameters(**args)

        # set the right background
        if self.transparency is False:
            if self.dark_mode:
                plt.style.use('dark_background')
            else:
                plt.style.use('default')
        # set the figure size
        self.gs = gridspec.GridSpec(2, 2)
        fig = plt.figure(figsize=self.fig_size)
        # choose the font
        if self.tex_font:
            if self.use_serif:
                # Serif latex font
                plt.rcParams.update({
                    "text.usetex": True,
                    "font.family": "serif",
                    "font.serif": ["Palatino"],
                }) 
            else:
                # Non-serif latex font
                plt.rcParams.update({
                    "text.usetex": True,
                    "font.family": "sans-serif",
                    "font.serif": ["Open Sans"],
                    "text.latex.preamble" : r"\usepackage{cmbright}"
                })
        
        self.fig = fig

    def add_panel(self, **args):
        self.update_parameters(**args)
        try:
            self.id_panel += 1
            if self.panel_position is None:
                self.ax.append(plt.subplot(self.n_line, self.n_colone, self.id_panel))
            else:
                i, j = self.panel_position
                if i is None:
                    self.ax.append(plt.subplot(self.gs[:, j]))
                elif j is None:
                    self.ax.append(plt.subplot(self.gs[i, :]))
                else:
                    self.ax.append(plt.subplot(self.gs[i, j]))
        except:
            self.ax = []
            self.id_panel = 1
            if self.panel_position is None:
                self.ax.append(plt.subplot(self.n_line, self.n_colone, self.id_panel))
            else:
                i, j = self.panel_position
                if i is None:
                    self.ax.append(plt.subplot(self.gs[:, j]))
                elif j is None:
                    self.ax.append(plt.subplot(self.gs[i, :]))
                else:
                    self.ax.append(plt.subplot(self.gs[i, j]))
        self.cpt_colors = 0

    def add_plot(self, **args):

        self.update_parameters(**args)

        # Pick color automatically
        if self.data_color is None:
            data_color = colorserie1[self.cpt_colors]
        elif isinstance(self.data_color, int):
            data_color = colorserie1[self.data_color]
        elif self.data_color == "autogray":
            if self.dark_mode:
                 data_color = np.array([0.9, 0.9, 0.9])
            else:
                data_color = np.array([0.1, 0.1, 0.1])
        else:
            data_color = self.data_color

        # In case of open symbol choice
        if self.open_symbols:
            self.markerfacecolor = 'none'
            if self.markeredgewidth == 0:
                self.markeredgewidth = 3
        else:
            self.markerfacecolor = data_color
            
        #assert self.x is not None
        self.ax[-1].plot(self.x,
                    self.y,
                    self.marker,
                    color = data_color,
                    markersize = self.markersize,
                    linewidth = self.data_linewidth,
                    label = self.data_label,
                    markeredgewidth = self.markeredgewidth,
                    markeredgecolor = data_color,
                    markerfacecolor = self.markerfacecolor)
        self.cpt_colors += 1

        # Convert the axis to log        
        if (self.type == 'semilogy') | (self.type == 'loglog'):
            self.ax[-1].set_yscale('log')
        if (self.type == 'semilogx') | (self.type == 'loglog'):
            self.ax[-1].set_xscale('log')

        self.reset_parameters()

    def add_subplotlabels(self, **args):
        """Add a labels to each axis of a figure.
        Four options:
        - For small letters, use type_label_panel='a'
        - For large letters, use type_label_panel='A'
        - Deactivate using type_label_panel=None
        - Provide a list of special character instead. The list must be the size of 'ax' 
        """
        self.update_parameters(**args)
        if self.type_label_panel is not None:
            if self.type_label_panel == "a":
                labels = []
                for i, value in zip(range(len(self.ax)),
                                    list(map(chr, range(97, 123)))):
                    if self.tex_font:
                        labels.append(r"$\textrm{"+str(value)+"}$")
                    else:
                        labels.append(value)
            elif self.type_label_panel == "A":
                labels = []
                for i, value in zip(range(len(self.ax)),
                                    list(map(chr, range(65, 91)))):
                    if self.tex_font:
                        labels.append(r"$\textrm{"+str(value)+"}$")
                    else:
                        labels.append(value)
            elif len(self.type_label_panel) >= 1:
                if len(self.type_label_panel) == len(self.ax):
                    labels = []
                    for value in self.type_label_panel:
                        labels.append(value)
                else: 
                    print("Inconsistent type_label_panel provided." \
                            "Parameter Ignored.")
                    labels = None
            else:
                labels = None
        else:
            labels = None

        if labels is not None:
            for i, subplotlabel in enumerate(labels):
                if self.sshift_label_panel is None:
                    trans = mtransforms.ScaledTranslation(
                        self.shift_label_panel, -0.2,
                        self.fig.dpi_scale_trans)
                else:
                    trans = mtransforms.ScaledTranslation(
                        self.sshift_label_panel[i], 0,
                        self.fig.dpi_scale_trans) 
                self.ax[i].text(
                    0.0,
                    1.0,
                    subplotlabel,
                    transform=self.ax[i].transAxes + trans,
                    va="top",
                    #bbox=dict(facecolor="none", alpha=0.5, edgecolor="none", pad=3.0),
                    fontdict = self.font,
                    color=self.axis_color
                )     


    def complete_panel(self, **args):
        
        self.update_parameters(**args)
        # Write label along x
        if self.xlabel is not None:
            self.ax[-1].set_xlabel(self.xlabel, fontdict=self.font)
            if self.cancel_x:
                self.ax[-1].set_xticklabels([])
        else:
            self.ax[-1].set_xticklabels([])
        # Write label along y
        if self.ylabel is not None:
            self.ax[-1].set_ylabel(self.ylabel, fontdict=self.font)
            if self.cancel_y:
                self.ax[-1].set_yticklabels([])  
        else:
            self.ax[-1].set_yticklabels([])
        # Write panel title
        if self.panel_title is not None:
            self.ax[-1].set_title(self.panel_title,
                                  fontdict=self.font,
                                  color=self.axis_color,
                                  y=1.02)

        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)
        self.ax[-1].yaxis.offsetText.set_fontsize(self.fontsize)
        self.ax[-1].minorticks_on()

        # Place tick (only in is implemented, and length is enforced)
        # Could be changed in the futur
        self.ax[-1].tick_params('both', length=10,
                                width=self.tickwidth1,
                                which='major', direction='in')
        self.ax[-1].tick_params('both', length=6,
                                width=self.tickwidth2,
                                which='minor', direction='in')
        self.ax[-1].xaxis.set_ticks_position('both')
        self.ax[-1].yaxis.set_ticks_position('both')

        # Apply linewidth to graph border
        for border in ["top", "bottom", "left", "right"]:
            self.ax[-1].spines[border].set_linewidth(self.panel_linewidth)

        # Estimate the right number of minor tick
        if self.locator_x == 'auto':
            scale = self.ax[-1].xaxis.get_scale()
            if scale == 'linear':
                locator_x = 2
            else:
                locator_x = None
        else:
            locator_x = self.locator_x
        #elif 
        if self.locator_y == 'auto':
            scale = self.ax[-1].yaxis.get_scale()
            if scale == 'linear':
                locator_y = 2
            else:
                locator_y = None   
        else:
            locator_y = self.locator_y   
        # Apply minor tick position
        if locator_x is not None:
            minor_locator_x = AutoMinorLocator(locator_x)
            self.ax[-1].xaxis.set_minor_locator(minor_locator_x)
        if locator_y is not None:
            minor_locator_y = AutoMinorLocator(locator_y)
            self.ax[-1].yaxis.set_minor_locator(minor_locator_y)

        if self.legend:
            if not self.ax[-1].get_legend_handles_labels() == ([], []):
                labels_in_fig = True
            else:
                labels_in_fig = False
            if labels_in_fig:
                self.ax[-1].legend(frameon=False, fontsize=self.fontsize,
                                   labelcolor=self.axis_color, loc='best',
                                   handletextpad=0.5, ncol=self.ncol_legend,
                                   handlelength = 0.86, borderpad = 0.3, 
                                   labelspacing=0.3)
                    
        # color the axis if requested
        if self.axis_color is not None:
            self.ax[-1].xaxis.label.set_color(self.axis_color)
            self.ax[-1].yaxis.label.set_color(self.axis_color)
            for axis in ["x", "y"]:
                self.ax[-1].tick_params(axis=axis, colors=self.axis_color)
                self.ax[-1].tick_params(axis=axis, which='both', colors=self.axis_color)
            for border in ["top", "bottom", "left", "right"]:
                self.ax[-1].spines[border].set_color(self.axis_color)
    
        if self.xpad is not None:
            self.ax[-1].tick_params(axis='x',
                                    colors=self.axis_color,
                                    pad = self.xpad)
        if self.ypad is not None:
            self.ax[-1].tick_params(axis='y',
                                    colors=self.axis_color,
                                    pad = self.ypad)

    def set_boundaries(self, **args):
        self.update_parameters(**args)
        if self.x_boundaries is not None:
            plt.xlim(self.x_boundaries)
        if self.x_ticks is not None:
            plt.xticks(self.x_ticks)
        if self.y_boundaries is not None:
            plt.ylim(self.y_boundaries)
        if self.y_ticks is not None:
            plt.yticks(self.y_ticks)  

    def save_figure(self, **args):
        self.update_parameters(**args)

        if self.saving_path is not None:
            assert os.path.exists(self.saving_path)
            saving_path = self.saving_path
            if saving_path[-1] != "/":
                saving_path = saving_path + "/"
        else:
            saving_path = "./"

        self.fig.tight_layout()
        if self.transparency is False:
            plt.style.use('default')
        if self.dark_mode:
            plt.savefig(saving_path + self.filename + "-dm.png",
                        bbox_inches = 'tight', pad_inches = 0.062,
                        transparent=self.transparency, dpi=200)
        else:
            plt.savefig(saving_path + self.filename + ".png",
                        bbox_inches = 'tight', pad_inches = 0.062,
                        transparent=self.transparency, dpi=200)
