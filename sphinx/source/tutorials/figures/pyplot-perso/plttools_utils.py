import os
import numpy as np

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

def random_lin_generator(xmin=0, xmax=10, slope=0.3, alpha=0.1, pref = 1, N = 50):
    x = np.linspace(xmin, xmax, N)
    y = slope*x + pref * (np.random.random(len(x))-0.5)*x**alpha
    return x, y

def multivariate_normal_distribution(meanx, meany, s1 = [0, 0], s2= [0, 0], N = 50):
    """Call the multivariate normal distribution of NumPy"""
    cov = np.array([s1, s2])
    pts = np.random.multivariate_normal([meanx, meany], cov, size=N)
    return pts[:,0], pts[:,1]

def random_log_generator(xmin=0, xmax=10, slope=0.3, alpha=0.1, pref = 1, N = 50):
    x = np.logspace(xmin, xmax, N)
    y = slope*x + pref * (np.random.random(len(x))-0.5)*x**alpha
    return x, y

def mygradient(N, color1, color2, final_value=False):
    """Generate a color gradient from color1 to color2"""
    if final_value:
        R = np.linspace(color1[0], color2[0], N)
        G = np.linspace(color1[1], color2[1], N)
        B = np.linspace(color1[2], color2[2], N)
    else:
        R = np.linspace(color1[0], color2[0], N)[1:-1]
        G = np.linspace(color1[1], color2[1], N)[1:-1]
        B = np.linspace(color1[2], color2[2], N)[1:-1]
    return np.vstack([R, G, B]).T
