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
from matplotlib import style
from matplotlib.patches import Ellipse
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
from hrvanalysis import plot_poincare
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - [%(levelname)-8s] (%(threadName)-10s) %(message)s")
PdfPages = PdfPages('results/monikit_results'+'.pdf')

def read_data(file,window=False):
    """
    reading the data files for RR intervals
    """
    data = pd.read_csv(file, delimiter=",", header=0)
    x = data["time_in_ms"] 
    #Define a timewindow based on the period
    if window:
        x2=x[0]+window*1000
        filtered_data = data[data.time_in_ms < x2 ]
        x = filtered_data["time_in_ms"]
        y = filtered_data["rr_in_ms"]
    else:
        y = data["rr_in_ms"]
    return x,y

def timeS2Poincare(rr_intervals,ndelay=1):
    ax1 = rr_intervals[:-ndelay]
    ax2 = rr_intervals[ndelay:]
    return ax1,ax2


def get_poincare_plot_features(rr_intervals):
    """
    Fuction to calculate
    SD1 is the standard deviation of projection of the Poincaré plot on the line perpendicular to the line of identity.
    SD2 is the standard deviation of the projection of the Poincaré plot on the line of identity (y=x).
    csi Cardiac Sympathetic Index
    cvi Cadiac Vagal Index
    
    Arguments
    ---------
    rr_intervals: list
    
    """
    diff_rr_intervals = np.diff(rr_intervals)
    # measures the std of transverse direction
    sd1 = np.sqrt(np.std(diff_rr_intervals, ddof=1) ** 2 * 0.5)
    # measures  std of longitudinal direction
    sd2 = np.sqrt(2 * np.std(rr_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_rr_intervals, ddof=1) ** 2)

    # measures the center of the poincare cloud
    mean_rri = np.mean(rr_intervals) 
    std_rri= np.std(rr_intervals) 
    
    # Measures the width and length of poincare cloud
    csi = sd2 / sd1
    T = 4 * sd1
    L = 4 * sd2
    cvi = np.log10(L * T)
    return sd1, sd2,csi, mean_rri,std_rri,cvi

def plot_timeseries_data(x,y,path=None,xlabel=None,ylabel=None,output=None,title="Rr Interval time series"):
    """
    The function plots the timeseries
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig()
    
def plot_distributions(rr_intervals,output=None,path=None,bin_size=100,title="Distribution of RR Intervals"):
    """
    Function plotting histogram distribution of the NN Intervals. Useful for geometrical features.
    
    Arguments
    ---------
    rr_intervals: list
        list of Normal to Normal Interval.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(title, fontsize=12)
    plt.xlabel("Time [ms]")
    plt.ylabel("Distribution of RR interval")
    plt.hist(rr_intervals,bin_size)
    #Add some statistics 
    sd1, sd2,csi, mean_rri,std,cvi= get_poincare_plot_features(rr_intervals)
    ax.text(0.90, 0.90, "mean= %.3f\n std=%.3f" % (mean_rri,std),
        horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))             
    plt.tight_layout()
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig() 
    
    
def plot_Poincare(rr_intervals,xlabel=None,ylabel=None,title="Poincaré Plot",output=None,path=None):
    """
    Fuction to plot Poincare plot with its features.
    
    Arguments
    ---------
    rr_intervals:List
    """
    x,y = timeS2Poincare(rr_intervals)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x,y,c='r', s=2)
    plt.xlim(min(rr_intervals) - 10, max(rr_intervals) + 10)
    plt.ylim(min(rr_intervals) - 10, max(rr_intervals) + 10)
    plt.title(title,fontsize=12)
    #plot the line IK=IK+1
    ax.plot(x,x,'black','-',markersize=10)
    
    sd1, sd2,csi, mean_rri,_,_= get_poincare_plot_features(rr_intervals)
    # Arrow plot settings
    sd1_arrow = ax.arrow(mean_rri, mean_rri, -sd1 * np.sqrt(2) / 2, sd1 * np.sqrt(2) / 2,
                         linewidth=3, ec='b', fc="b", label="SD1")
    sd2_arrow = ax.arrow(mean_rri, mean_rri, sd2 * np.sqrt(2) / 2, sd2 * np.sqrt(2) / 2,
                         linewidth=3, ec='g', fc="g", label="SD2")
    plt.legend(handles=[sd1_arrow, sd2_arrow], fontsize=12, loc="best")
    ax.text(0.80, 0.90, "CSI= %.3f" % (csi),
        horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig()


"""
Question 1
"""
#loading  data from the file data.csv
x,rr_intervals_list =read_data(file='data/data.csv')
#plotting data
plot_timeseries_data(x,rr_intervals_list,xlabel='time [ms]',ylabel='RR interval[ms]',
                     output='1.data', path='results/')

"""
Question 2
"""

#Calculate statistics(a)
plot_distributions(rr_intervals_list,
                   output='2.a_distribution', path='results/')
#plotting showPoincare plot with csi (b)
plot_Poincare(rr_intervals_list,xlabel='RRn[ms]',ylabel='RRn+1[ms]',
              output='2.b_Poincare',path='results/')
"""
Question 3
"""
#loading  data from the file data.csv with a window size of 30 seconds
x2,rr_intervals_list2 =read_data(file='data/data.csv',window=30)
#plotting data of 30 seconds window
plot_timeseries_data(x2,rr_intervals_list2,xlabel='time [ms]',ylabel='RR interval[ms]',
                     output='3.data_30sec', path='results/',
                     title="Rr Interval time series in 30 sec window")
#Calculate statistics(a)
plot_distributions(rr_intervals_list2,title="Distribution of RR Intervals (Window size of 30 sec)",
                   output='3.a_distribution_30sec', path='results/',bin_size=10)
#plotting showPoincare plot with csi (b)
plot_Poincare(rr_intervals_list2,xlabel='RRn[ms]',ylabel='RRn+1[ms]',title="Poincaré Plot (Window size of 30 sec)",
              output='3.b_Poincare_30sec',path='results/')
PdfPages.close()
