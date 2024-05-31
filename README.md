# ECG Signal Analysis for Seizure Detection

## Introduction
This script focuses on analyzing RR intervals derived from ECG signals to detect epileptic seizures. RR intervals represent the time intervals between successive R-peaks, which correspond to individual heartbeats. By studying variations in RR intervals, we aim to identify patterns associated with epileptic seizures.

## Project Objective
The primary objective of this project is to develop algorithms for automated seizure detection based on ECG signal analysis. By detecting specific patterns in RR interval variations indicative of seizure events, we aim to create reliable tools for seizure prediction and monitoring.

## Approach
The approach involves the following steps:
1. **Data Collection**: Obtain ECG data containing RR intervals from both healthy individuals and epileptic patients.
2. **Feature Extraction**: Implement statistical and signal processing techniques to extract relevant features from RR interval data.
3. **Pattern Recognition**: Analyze RR interval variations to identify distinct patterns associated with epileptic seizures.
4. **Algorithm Development**: Develop algorithms for automated seizure detection based on identified patterns.
5. **Evaluation**: Evaluate the performance of the developed algorithms using validation datasets and clinical studies.

## Features
- **Mean RR Interval**: Average duration between consecutive R-peaks.
- **Standard Deviation of RR Intervals**: Measure of variability in RR interval durations.
- **Cardiosympathetic Index (CSI)**: Derived metric indicating cardiac sympathetic activity.

## Data Description
- The dataset contains ECG data in the form of RR intervals.
- Seizure events are annotated in the dataset, including the type of seizure and the onset time.
- Additional metadata provides information about the recording source and patient demographics.

## Installation and Usage
1. **Clone Repository**: Clone the repository to obtain the source code.
```
git clone git@github.com:ahmedqamesh/monikit_assignment.git
```

2. **Dependencies**: Install required dependencies using pip.
```
pip install hrvanalysis
```
3. **Run the Code**: Execute the main script to perform ECG signal analysis and seizure detection.
```
python main_assignment_test.py
```

### Description of the data file:
the data file data.zip contains a scientific paper and three comma-separated text files:
- data.csv contains two columns: the timestamp of each RR interval and the length of this interval. Both given in milliseconds. This data is from one continuous record of about 160 hours.
- data_header.csv contains general information about the record which was the source for the data
- seizures.csv contains information about the seizures occurred in the given record.
Important columns are:
- classification which can be ‘SP’ for simple partial, ‘CP’ for complex partial and ‘UC’ for unclassified
- eeg_onset_secs which gives the beginning of the seizure given in seconds from start of the record All other columns can be ignored. Refer to the seizures only by IDs given by their row number (1:22) and ignore the seizure_id column.


## Contributing and Contact Information:
We welcome contributions from the community please contact : `ahmed.qamesh@gmail.com`.
