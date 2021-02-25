clear, clc, close all

% get a section of the sound file
[x, fs] = audioread('track.wav');   % load an audio file
x = x(:, 1);                        % get the first channel
x = x/max(abs(x));                  % normalize the signal
N = length(x);                      % signal length
t = (0:N-1)/fs;                     % time vector

% plot the signal
figure(1)
hPlot = plot(t, x, 'r');
xlim([0 max(t)])
ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Normalized amplitude')
title('The signal in the time domain')

% enable interactive section selection
secsel(hPlot)

% Instructions for use:
% - Double click the left mouse button to place the first mark line to a 
%   new possition;
% - Double click the right mouse button to place the second mark line to a 
%   new possition;
% - Press "Enter" key in order to save the selected 2D plot section as a
%   mat. data file;
% - Press "End" key in order to close the interactive section selection
%   regime.

