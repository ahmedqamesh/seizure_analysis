### Package task:
1. Load the data and plot it. What do you generally notice when you look at the data?
Describe it shortly.
2. Implement the following three features:
a. Mean
b. Standard deviation
c. Cardiosympathetic Index (CSI) as described in the given paper
3. Calculate the features on the given data with a window size of 30 seconds and overlap
between windows of 5 seconds
4. Make plots of the features for all seizures where the data is good enough
5. Do you think these features are helpful for our main goal, the detection of epileptic
seizures? What would you do next?
6. Prepare a 15-20 Minute presentation with all your outcomes and present it to us during
the interview. Choose whatever Presentation-Tool you want (PowerPoint, Keynote, Pen
and Paper,...).

### Description of the data file:
the data file data.zip contains a scientific paper and three comma-separated text files:
- data.csv contains two columns: the timestamp of each RR interval and the length of this interval. Both given in milliseconds. This data is from one continuous record of about 160 hours.
- data_header.csv contains general information about the record which was the source for the data
- seizures.csv contains information about the seizures occurred in the given record.
Important columns are:
- classification which can be ‘SP’ for simple partial, ‘CP’ for complex partial and ‘UC’ for unclassified
- eeg_onset_secs which gives the beginning of the seizure given in seconds from start of the record All other columns can be ignored. Refer to the seizures only by IDs given by their row number (1:22) and ignore the seizure_id column.


### Installation and usage

Clone the repository to get a copy of the source code (for developers):

```
git clone git@github.com:ahmedqamesh/monikit_assignment.git
```

In order to run the package use the following command
```
python main_assignment_test.py
```

### Dependencies:
- hrvanalysis
- 
