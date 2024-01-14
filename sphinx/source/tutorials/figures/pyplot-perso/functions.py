import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.ticker import AutoMinorLocator

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})

fontsize = 34
font = {'family': 'sans', 'color':  'black', 'weight': 'normal', 'size': fontsize}


def prepare_figure(mode, transparency = False, desired_figsize=(18,6)):
    if transparency is False:
        if mode == 'dark':
            plt.style.use('dark_background')
        elif mode == 'light':
            plt.style.use('default')
        else:
            print("WARNING: unknown choice of mode")
    fig = plt.figure(figsize=desired_figsize)
    return fig

def add_subplotlabels(fig, ax, labels, shift=0.2, specific_shift=None, color=None):
    """Add a labels to each axis of a figure."""
    assert len(ax) == len(labels), """WARNING: number of labels different from the number of subplots"""

    for i, subplotlabel in enumerate(labels):
        if specific_shift is None:
            trans = mtransforms.ScaledTranslation(
                shift, -0.2, fig.dpi_scale_trans)
        else:
            trans = mtransforms.ScaledTranslation(
                specific_shift[i], 0, fig.dpi_scale_trans) 
        if color is None:
            ax[i].text(
                0.0,
                1.0,
                subplotlabel,
                transform=ax[i].transAxes + trans,
                va="top",
                #bbox=dict(facecolor="white", alpha=0.5, edgecolor="none", pad=3.0),
                fontdict = font,
            )
        else:
            ax[i].text(
                0.0,
                1.0,
                subplotlabel,
                transform=ax[i].transAxes + trans,
                va="top",
                #bbox=dict(facecolor="none", alpha=0.5, edgecolor="none", pad=3.0),
                fontdict = font,
                color=color
            )     


def complete_panel(ax, xlabel, ylabel, cancel_x=False, cancel_y=False,
                   font=font, fontsize=fontsize, linewidth=2.5, tickwidth1=2.5,
                   tickwidth2=2, legend=True, ncol=1, locator_x = 2, locator_y = 2,
                   title=None, axis_color=None, xpad = None, ypad = None):
    
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontdict=font)
        if cancel_x:
            ax.set_xticklabels([])
    else:
        ax.set_xticklabels([])

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontdict=font)
        if cancel_y:
            ax.set_yticklabels([])  
    else:
        ax.set_yticklabels([])

    if title is not None:
        ax.set_title(title, fontdict=font)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.minorticks_on()

    ax.tick_params('both', length=10, width=tickwidth1, which='major', direction='in')
    ax.tick_params('both', length=6, width=tickwidth2, which='minor', direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    # border of graph
    ax.spines["top"].set_linewidth(linewidth)
    ax.spines["bottom"].set_linewidth(linewidth)
    ax.spines["left"].set_linewidth(linewidth)
    ax.spines["right"].set_linewidth(linewidth)

    minor_locator_x = AutoMinorLocator(locator_x)
    ax.xaxis.set_minor_locator(minor_locator_x)
    minor_locator_y = AutoMinorLocator(locator_y)
    ax.yaxis.set_minor_locator(minor_locator_y)

    if legend:
        ax.legend(frameon=False, fontsize=fontsize, labelcolor=axis_color,
                loc='best', handletextpad=0.5, ncol=ncol,
                handlelength = 0.86, borderpad = 0.3, 
                labelspacing=0.3)
                
    if axis_color is not None:
        ax.xaxis.label.set_color(axis_color)
        ax.yaxis.label.set_color(axis_color)
        ax.tick_params(axis='x', colors=axis_color)
        ax.tick_params(axis='y', colors=axis_color)
        ax.spines['left'].set_color(axis_color)
        ax.spines['top'].set_color(axis_color)
        ax.spines['bottom'].set_color(axis_color)
        ax.spines['right'].set_color(axis_color)
        ax.tick_params(axis='y', which='both', colors=axis_color)
        ax.tick_params(axis='x', which='both', colors=axis_color)

    if xpad is not None:
        ax.tick_params(axis='x', colors=axis_color, pad = xpad)

    if ypad is not None:
        ax.tick_params(axis='y', colors=axis_color, pad = ypad)

def save_figure(plt, fig, mode, git_root, path_figures, filename, show=False, transparency=True):
    assert os.path.exists(git_root + path_figures)
    fig.tight_layout()
    if mode == 'light':
        if transparency is False:
            plt.style.use('default')
        plt.savefig(git_root + path_figures + filename + "-light.png",
                    bbox_inches = 'tight', pad_inches = 0.062,
                    transparent=transparency, dpi=200)
    else:
 
        plt.savefig(git_root + path_figures + filename + "-dark.png",
                    bbox_inches = 'tight', pad_inches = 0.062,
                    transparent=transparency, dpi=200)
    if show:
        plt.show()

def set_boundaries(plt, x_boundaries=None, x_ticks=None, y_boundaries=None, y_ticks=None):
    if x_boundaries is not None:
        plt.xlim(x_boundaries)
    if x_ticks is not None:
        plt.xticks(x_ticks)
    if y_boundaries is not None:
        plt.ylim(y_boundaries)
    if y_ticks is not None:
        plt.yticks(y_ticks)  

def import_ave_time(filename, folder=None):
    assert filename[:6] == "output"
    if folder is None:
        if os.path.exists(filename):
            data = np.loadtxt(filename)
            try:
                time, data = data.T
            except:
                time, data, _ = data.T
            if os.path.exists("data_plot/") is False:
                os.mkdir("data_plot/")
            np.savetxt("data_plot/"+filename[7:], np.vstack([time, data]).T)
        else:
            time, data = np.loadtxt("data_plot/"+filename[7:]).T
    else:
        if os.path.exists(folder+filename):
            data = np.loadtxt(folder+filename)
            try:
                time, data = data.T
            except:
                time, data, _ = data.T
            if os.path.exists("data_plot/") is False:
                os.mkdir("data_plot/")
            np.savetxt("data_plot/"+filename[7:], np.vstack([time, data]).T)
    return time, data