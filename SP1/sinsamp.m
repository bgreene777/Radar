% Brian R. Greene
% sinsamp.m
%
% This program will sample a sinusoidal signal
% and plot the result.
%
clc
clear
close all

% freq = 20 Hz, samp = 100 Hz
samptime=1/100; % s
freq=20; % Hz
time=[0:99]*samptime;
x=cos(2*pi*freq*time);
subplot(2,1,1);
plot(time,x);
title('Sinusoid sampled at 100 Hz');
xlabel('Time (sec)');
ylabel('x(t)');
grid on;

%
% freq = 20 Hz, samp = 20 Hz
%
clear all;
samptime=1/20;
freq=20;
time=[0:19]*samptime;
x=cos(2*pi*freq*time);
subplot(2,1,2);
plot(time,x);
title('Sinusoid sampled at 20 Hz');
xlabel('Time (sec)');
ylabel('x(t)');
grid on;