function plotData(x, y)
%PLOTDATA Plots the data points x and y into a new figure 
%   PLOTDATA(x,y) plots the data points and gives the figure axes labels of
%   the timestamp of each RR interval and the length of this interval.

fig = figure(); % open a new figure window
plot(x, y, 'rx', 'MarkerSize', 10); % Plot the data
ylabel('RR interval[ms]'); % Set the yaxis label
xlabel('time [ms]'); % Set the xaxis label
title ("RR data of 160 hours");
print (fig, "data.pdf");
% ============================================================
end
