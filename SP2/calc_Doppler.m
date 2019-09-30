% Brian R. Greene
% Plot I(t) and Q(t) together for several aximuth angles and range gates

clc
clear
close all

%% Load data
load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200821_e02.60.mat');

%% Get actual r values
r_all = linspace(0, 1009, 1010) * delr / 1000;

%% Inside storm: range gate 600, azimuth angle 2
r_1 = r_all(600);
az_1 = az_set(2);

% X_h
f1=figure(1);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_h(2,600,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_h(2,600,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_h at az = ', num2str(az_1, '%3.0f'), '^\circ, r = ', num2str(r_1, '%6.3f'), ' km'])
% save
% print(f1, 'Q6_storm_h', '-dpng', '-r300')