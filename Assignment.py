#
# ------------------------------------------------------------
# Copyright (c) All rights reserved
# Ahmed Qamesh, Institute of Physics, University of Bonn
# 15/05/2019 14:50
# ------------------------------------------------------------
#
'''
    Script to solve Assignment  “First insights in the data”
    for the Full-Time-Position Data Scientist
    at monikit
    For more information, check: monikit Assignment Data Scientist.pdf
'''
from __future__ import division
from numpy import loadtxt, arange, append
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import csv
import tables as tb
from matplotlib import style
from matplotlib.patches import Ellipse
from matplotlib import pyplot as p
from mpl_toolkits.mplot3d import Axes3D    # @UnusedImport
from math import pi, cos, sin
import logging
import os
import time
import matplotlib as mpl
from matplotlib import colors
import pylab as P
import pandas as pd
from hrvanalysis import plot_poincare
from matplotlib.animation import FuncAnimation
from matplotlib import animation, rc
from matplotlib.pyplot import xticks
from matplotlib import gridspec
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - [%(levelname)-8s] (%(threadName)-10s) %(message)s")
logger = logging.getLogger(__name__)
PdfPages = PdfPages('results/monikit_results'+'.pdf')

def read_data(data_file):
    """
    reading the data files for RR intervals
    Argument
    ---------
    data_file: file to be read
    """
    data = pd.read_csv(data_file, delimiter=",", header=0)
    x = data["time_in_ms"] 
    y = data["rr_in_ms"]
    logger.info('Loading data file')
    return x,y

def read_seizures(seizures_file):
    """
    reading the seizures file
    Argument
    ---------
    seizures_file: file to be read
    """
    seizures = pd.read_csv(seizures_file, delimiter=",", header=0)
    classification = seizures["classification"]
    eeg_onset_secs= seizures["eeg_onset_secs"]
    logger.info('Loading seizures file')
    return classification,eeg_onset_secs
    
def timeS2Poincare(rr_intervals,ndelay=1):
    """
    create a delay of ndelay in rr_intervals 
    Argument
    --------- 
    rr_intervals: list 
    """
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
    
    References
    ----------
    .. [5] Using Lorenz plot and Cardiac Sympathetic Index of heart rate variability for detecting \
    seizures for patients with epilepsy, Jesper Jeppesen et al, 2014
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

def get_rolling_window_features(file,time_winow=30,overlap=5):
    """
    reading the data files for RR intervals and calculate the features based on the overlapping parameters (time_winow and ovelap).
    Arguments
    ---------
    file:list
    time_winow: int
    overlap : int
    """
    logger.info('Calculating the features on the data within a time window of %i s and ovelapping of %i s'%(time_winow,overlap))
    # Reshape the data with window size and overlapping.
    data = pd.read_csv(file, delimiter=",", header=0)
    overlapped_data=windows_overlap(data,time_winow*1000,overlap*1000)
    rolled_window_mean=overlapped_data.rolling(time_winow*1000).mean()
    rolled_window_std=overlapped_data.rolling(time_winow*1000).std()
    #plot the data after window
    return rolled_window_mean, rolled_window_std, overlapped_data
    
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


    
def get_csi_overlap(file,time_winow=30,overlap=5):
    """
    The function returns csi for the overlapped data
    Arguments
    ---------
    file:list
    time_winow: int
    overlap : int
    """
    logger.info('Calculating CSI on the data within a time window of %i s and ovelapping of %i s'%(time_winow,overlap))
    data = pd.read_csv(file, delimiter=",", header=0)
    overlapped_data=windows_overlap(data,time_winow*1000,overlap*1000)
    #csi_overlapped = overlapped_data["rr_in_ms"].rolling(time_winow*1000).apply(get_features,raw=True,args="rr_in_ms")
    rr_in_ms = overlapped_data["rr_in_ms"]
    csi_overlapped=[]
    time_in_chunk_ms = []
    def get_csi_feature(rr_intervals):
        """
        The function returns csi for the overlapped data
        Arguments
        ---------
        rr_intervals:list
        """    
        diff_rr_intervals = np.diff(rr_intervals)
        # measures the std of transverse direction
        sd1 = np.sqrt(np.std(diff_rr_intervals, ddof=1) ** 2 * 0.5)
        # measures  std of longitudinal direction
        sd2 = np.sqrt(2 * np.std(rr_intervals, ddof=1) ** 2 - 0.5 * np.std(diff_rr_intervals, ddof=1) ** 2)
        csi = sd2 / sd1
        return csi
    
    for i in range(0,len(rr_in_ms),time_winow*1000):
        rr_in_ms_chunk =rr_in_ms[i:i+time_winow*1000]
        time_in_chunk_ms =np.append(time_in_chunk_ms,i+time_winow)  
        csi_overlapped =np.append(csi_overlapped,get_csi_feature(rr_in_ms_chunk))
    return csi_overlapped, time_in_chunk_ms

def analyse_seizures(seizures_file="data/seizures.csv",data_file='data/data.csv',interval=10000):
    """
    The function to plot features for all seizures where the data is good enough,
    Argument
    ---------   
    seizures_file: list
    data_file: list
    interval: int
        a window in [ms] of data to cover a time line in ms for the data collected around the seizure time
    """
    #read the seizures file 
    classification,eeg_onset_secs = read_seizures(seizures_file)
    #index elemets in the seizures file 
    get_indexes = lambda eeg_onset_secs, xs: [i for (y, i) in zip(xs, range(len(xs))) if eeg_onset_secs == y]
    
    #filter the data based on the seizures
    data = pd.read_csv(data_file, delimiter=",", header=0)
    index = 0
    csi_array = []
    std_array = []
    mean_array = []
    classification_array=[]
    eeg_onset_array=[]
    for eeg in eeg_onset_secs:
        data_filtered= data[(data.time_in_ms <= eeg*1000+interval) & (data.time_in_ms >=eeg*1000-interval)]
        x = data_filtered["time_in_ms"] 
        rr_intervals_list = data_filtered["rr_in_ms"]
        
        #plot each seizure 
        seizure_classification = classification[index]
        _, _,csi, mean_rri,std_rri,_ = get_poincare_plot_features(rr_intervals_list)
        csi_array = np.append(csi_array, csi)
        std_array = np.append(std_array, std_rri)
        mean_array = np.append(mean_array, mean_rri)
        classification_array = np.append(classification_array, seizure_classification)
        eeg_onset_array = np.append(eeg_onset_array, eeg_onset_secs[index])
        plot_timeseries_data_seizures(x, rr_intervals_list, path='results/', xlabel="Time [ms]",
                                       ylabel="RR interval [ms]", output="4.data_"+str(index+1)+"seizure",
                                       title="Rr interval time series",scale=1000, classification= seizure_classification, 
                                       eeg_onset_secs=eeg*1000,index=index+1, csi=csi, mean_rri=mean_rri,std_rri=std_rri)
        index += 1
    return mean_array, std_array, csi_array,classification_array,eeg_onset_array


def plot_timeseries_data_seizures(x,y,path=None,xlabel=None,ylabel=None,output=None,title="Rr Interval time series",
                         add_seizures=False,scale=1000,classification=None,eeg_onset_secs=None,index=None,csi=None, mean_rri=None,std_rri=None):
    """
    The function plots the timeseries around each seizure
    """
    style.use("seaborn-darkgrid")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    plt.axvline(x=eeg_onset_secs, linewidth=2, color='#d62728', linestyle='dashed')
    label = "seizure:"+str(index)+" "+classification+"\n >>Features:\n mean=%0.3f\n std=%0.3f\n CSI=%0.3f"%(mean_rri,std_rri,csi)
    ax.text(0.95, 0.95, label,
            horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path+ output +".png", bbox_inches='tight')
    PdfPages.savefig()
    
def plot_timeseries_data(x,y,path=None,xlabel=None,ylabel=None,output=None,title="Rr Interval time series",
                         add_seizures=False,scale=1000,classification=False,):
    """
    The function plots the timeseries
    """
    logger.info('plotting the timeseries data')
    style.use("seaborn-darkgrid")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    if add_seizures:
        classification,eeg_onset_secs = read_seizures("data/seizures.csv")
        for i in eeg_onset_secs:
            plt.axvline(x=i*scale, linewidth=1, color='#d62728', linestyle='dashed')
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
    logger.info('Plotting histogram distribution of the NN Intervals')
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
    logger.info('Creating Poincare plot with its features')
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
        
def plot_seizure_features(data1,data2,data3,output=None,path=None,classification=None,eeg_onset_array=None,interval=50000):
    """
   plot the data the features based on the overlapping parameters (time_winow and ovelap).
    """  
    logger.info('Plotting seizure features')
    fig =plt.figure()
    ax1 = fig.add_subplot(111)
    plt.title("data features within a window of %i ms around each feature"%interval,fontsize=12)
    mean_ax = ax1.plot(data1,label="mean")
    Std_ax = ax1.plot(data2,label="Std")
    CSI_ax = ax1.plot(data3,label="CSI")
    ax1.set_xlabel("seizures ID")
    plt.xticks(np.arange(0,len(data1)), classification)
    plt.legend(loc="upper right", prop={'size': 12})
    plt.tight_layout()
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig() 

def plot_rolling_window_features(data1,data2,xlabel=None,ylabel=None,
                                 time_winow=30,overlap=5,output=None,path=None):
    """
   plot the data the features based on the overlapping parameters (time_winow and ovelap).
    """  
    logger.info('Plotting features with a rolling window of %i s (overlapping %i s)'%(time_winow,overlap))
    ticks = np.arange(0,len(data1),len(data1)/8)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1 = data1["rr_in_ms"].plot(x="Time [ms]",legend=True,label="Mean",rot=45,title ="Data averages with window size of %i s (overlapping %i s)"%(time_winow,overlap))
    ax2 = data2["rr_in_ms"].plot(secondary_y=True,legend=True,label="Std",rot=45)
    ax.set_ylabel("Mean")
    mean_mean=np.mean(data1["rr_in_ms"])
    mean_std=np.mean(data2["rr_in_ms"])
    ax.text(0.77, 0.80, "mean= %.3f\n std=%.3f" % (mean_mean,mean_std),
        horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7)) 
    ax.right_ax.set_ylabel('Std')
    ax.set_xlabel("Time [ms]")
    plt.tight_layout()
    plt.savefig(path+ output + ".png", bbox_inches='tight')
    PdfPages.savefig()
        
"""
   Question 1
"""
logger.info('######################[Question 1]#####################################')
   #loading  data from the file data.csv
x,rr_intervals_list =read_data(data_file='data/data.csv') 
   #plotting  data from the file data.csv
plot_timeseries_data(x,rr_intervals_list,path='results/',xlabel="Time [ms]",ylabel="RR interval [ms]",output="1.data",title="Rr interval time series")

"""
Question 2
"""
logger.info('######################[Question 2]#####################################')
 #Calculate statistics(a)
plot_distributions(rr_intervals_list, output='2.a_distribution', path='results/')
#plotting showPoincare plot with csi (b)
plot_Poincare(rr_intervals_list,xlabel='RRn[ms]',ylabel='RRn+1[ms]',
                output='2.b_Poincare',path='results/')
"""
 Question 3
"""
logger.info('######################[Question 3]#####################################')
 #Calculating the features on the data within a time window of 30s and ovelapping of 5s
rolled_window_mean, rolled_window_std,overlapped_data = get_rolling_window_features(file='data/data.csv',time_winow=30,overlap=5)

#Calculating csi (b) for a time window of 30s and ovelapping of 5s
plot_rolling_window_features(rolled_window_mean, rolled_window_std,output='3.data_features_30sec', path='results/',time_winow=30,overlap=5) 

csi_overlapped, time_in_chunk_sec = get_csi_overlap('data/data.csv')
#Plotting csi (b) for a time window of 30s and ovelapping of 5s
plot_timeseries_data(time_in_chunk_sec,csi_overlapped,path='results/',xlabel="Time [ms]",add_seizures=False,scale=1,
                      ylabel="Cardiac Sympathetic Index(CSI)",output="3.data_csi_30sec",title="CSI for data with window size of 30s (overlapping 5s)")
"""
 Question 4
 """
logger.info('######################[Question 4]#####################################')
   #plotting  data and the corresponfing seizures
plot_timeseries_data(x,rr_intervals_list,path='results/',xlabel="Time [ms]",ylabel="RR interval [ms]",output="4.seizures_data",
                     title="Rr interval time series [The red lines represent the seizures]",add_seizures=True,scale=1000)
#Calculating the features on the data arround each seizure
mean_array, std_array, csi_array,classification, eeg_onset_array= analyse_seizures(seizures_file="data/seizures.csv",data_file='data/data.csv',interval=50000)
#plotting the features on the data arround each seizure
plot_seizure_features(mean_array, std_array, csi_array,output='4.seizure_features',path='results/',classification=classification,eeg_onset_array=eeg_onset_array,interval=50000)
PdfPages.close()