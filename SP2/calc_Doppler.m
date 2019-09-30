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

% plot vertical lines for beginning and end of 1 cycle in I
x1 = 0.003;
x2 = 0.0105;
y1=get(gca,'ylim');
plot([x1 x1], y1, 'k')
plot([x2 x2], y1, 'k')

% plot vertical dashed lines for beginning and end of 1 cycle in H
x3 = 0.01;
x4 = 0.0165;
y1=get(gca,'ylim');
plot([x3 x3], y1, '--k')
plot([x4 x4], y1, '--k')

% save
print(f1, 'Q7_IQ', '-dpng', '-r300')

%% 1 cycle I
T = x2 - x1;
f_d = 1/T;
% use f_d = -2 * v_r / lambda
v_r = -f_d * lambda / 2;

%% 1 cycle Q
T2 = x4 - x3;
f_d_2 = 1/T2;
v_r_2 = -f_d_2 * lambda / 2;