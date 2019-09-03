% Brian R. Greene
% complexsin.m
%
% This program will create a complex sinusoid
%
clc
clear
close all

t = 0:0.001:1;
len = length(t);
v=10*sin(2*pi*10*t) + j*10*cos(2*pi*10*t+0.3) +randn(1,len) + j*randn(1,len);

%% Plot
f1 = figure(1);
f1.Position = [375,492,1100,400];

plot(t, v)
ax1 = gca;
ax1.FontSize = 14;
title('Real part of complex vector vs. time', 'FontSize', 20)
xlabel('t [s]', 'FontSize', 20)
ylabel('Re(v)', 'FontSize', 20)
grid on

% save
print(f1, 'Q4', '-dpng', '-r200')