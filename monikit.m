%  Instructions:

%% Initialization
clear ; close all; clc
%% ======================= Part 1: Plotting =======================
fprintf('Plotting Data ...\n')
data = csvread('data/data.csv');
x = data(1:100, 1); y = data(1:100, 2);
m = length(y); % number of training examples
% Plot Data
plotData(x, y);
fprintf('Program paused. Press enter to continue.\n');
pause;
%% ======================= Part 1: Statistics =======================
mean =mean(y);
std = std(y);
