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
from matplotlib import animation
import pandas as pd
from hrvanalysis import plot_poincare
from matplotlib.animation import FuncAnimation
from matplotlib import animation, rc
from IPython.display import HTML, Image
from numpy.lib.stride_tricks import as_strided
from matplotlib.pyplot import xticks
# equivalent to rcParams['animation.html'] = 'html5'
rc('animation', html='html5')

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - [%(levelname)-8s] (%(threadName)-10s) %(message)s")
PdfPages = PdfPages('results/monikit_results'+'.pdf')

def read_data(file):
    """
    reading the data files for RR intervals
    """
    data = pd.read_csv(file, delimiter=",", header=0)
    x = data["time_in_ms"] 
    y = data["rr_in_ms"]
    return x,y

def read_seizures(file):
    """
    reading the seizures files
    """
    data = pd.read_csv(file, delimiter=",", header=0)
    classification = data["classification"]
    eeg_onset_secs= data["eeg_onset_secs"]
    print (classification,eeg_onset_secs)
    return classification,eeg_onset_secs

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
    style.use("seaborn-darkgrid")
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
    style.use("seaborn-darkgrid")
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
    style.use("seaborn-darkgrid")
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
        
def plot_rolling_window_features(data1,data2,xlabel=None,ylabel=None,title="Poincaré Plot",output=None,path=None):
    
    style.use("seaborn-darkgrid") 
    ticks = np.arange(0,len(data1),len(data1)/8)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1 = data1["rr_in_ms"].plot(x="Time [ms]",legend=True,label="Mean",rot=45,title ="Data with window size of 30s (overlapping 5s)")
    ax2 = data2["rr_in_ms"].plot(secondary_y=True,legend=True,label="Std",rot=45)
    ax.set_ylabel("Mean")
    mean_mean=np.mean(data1["rr_in_ms"])
    mean_std=np.mean(data2["rr_in_ms"])
    ax.text(0.80, 0.80, "mean= %.3f\n std=%.3f" % (mean_mean,mean_std),
        horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7)) 
    ax.right_ax.set_ylabel('Std')
    ax.set_xlabel("Time [ms]")
    plt.tight_layout()
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig() 
    
def windows_overlap(d, w, t):  
    """
    The function returns a new dataframe with overlapping
    Arguments
    ---------
    d: dataframe (list of data)
    w: window size(eg. 30 sec) 
    t: overlapping factor (eg. 5sec of window size 30 sec).
    """
    r = np.arange(len(d))   
    s = r[::t]   
    z = list(zip(s, s + w))   
    f = '{0[0]}:{0[1]}'.format
    g = lambda t: d.iloc[t[0]:t[1]]   
    return pd.concat(map(g, z), keys=map(f, z)) 

#plot_data_animation()
def get_rolling_window_features(file,time_winow=30,overlap=5):
    """
    reading the data files for RR intervals.
    Arguments
    ---------
    time:time series
    data: data collected for rr_intervals_list
    """
    # Reshape the data with window size and overlapping.
    data = pd.read_csv(file, delimiter=",", header=0)
    overlapped_data=windows_overlap(data,time_winow*1000,overlap*1000)
    rolled_window_mean=overlapped_data.rolling(time_winow*1000).mean()
    rolled_window_std=overlapped_data.rolling(time_winow*1000).std()
    # plot the data after window
    return rolled_window_mean, rolled_window_std, overlapped_data
  
# """
# Question 1
# """
# #loading  data from the file data.csv
# x,rr_intervals_list =read_data(file='data/data.csv') 
# """
# Question 2
# """
#  #Calculate statistics(a)
# plot_distributions(rr_intervals_list,
#                     output='2.a_distribution', path='results/')
#  #plotting showPoincare plot with csi (b)
# plot_Poincare(rr_intervals_list,xlabel='RRn[ms]',ylabel='RRn+1[ms]',
#                output='2.b_Poincare',path='results/')
# """
# Question 3
# """
# #loading  data from the file data.csv with a window size of 30 seconds
# rolled_window_mean, rolled_window_std,overlapped_data = get_rolling_window_features(file='data/data.csv')
# plot_rolling_window_features(rolled_window_mean, rolled_window_std,output='3.data_features_30sec', path='results/')
# 
# #plotting showPoincare plot with csi (b)
# plot_Poincare(overlapped_data["rr_in_ms"],xlabel='RRn[ms]',ylabel='RRn+1[ms]',title="Poincaré Plot (Window size of 30 sec)",
#               output='3.b_Poincare_30sec',path='results/')
# 
# PdfPages.close()
"""
Question 4
"""
read_seizures(file='data/seizures.csv') 

