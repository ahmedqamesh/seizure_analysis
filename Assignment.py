from __future__ import division
from numpy import loadtxt, arange
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mtick
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.backends.backend_pdf import PdfPages
import csv
from scipy.optimize import curve_fit
import tables as tb
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
import itertools
from matplotlib.colors import LogNorm
from matplotlib import pyplot as p
from mpl_toolkits.mplot3d import Axes3D    # @UnusedImport
from math import pi, cos, sin
import logging
from scipy.linalg import norm
import os
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib.transforms as mtransforms
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import pylab as P
import pandas as pd
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - [%(levelname)-8s] (%(threadName)-10s) %(message)s")

global PdfPages
def read_data1(path):
    with open(path+".csv", 'r')as data:
        x,y =[],[]
        reader = csv.reader(data)
        reader.next()
        for row in reader:
            x = np.append(x, np.float(row[0]))  # time_in_ms
            y = np.append(y, np.float(row[1]))  # RR_in_ms
    return x,y


def read_data(path):
    data = pd.read_csv(path)
    return data.values

def plot_data(x,y,path=None,xlabel=None,ylabel=None,output=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(x, y, fmt='o', markersize='4')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(path + ".png", bbox_inches='tight')
    plt.tight_layout()
    plt.show()
    PdfPages.savefig()
    
    
x,y =read_data('data/data')
#plot = plot_data (x,y, path='data/data',xlabel='time [ms]',ylabel='RR interval[ms]',output='data')
PdfPages = PdfPages('monikit_Bonn'+'.pdf')
PdfPages.close()