% Brian R. Greene
% METR 5673 HW2, DZ Q3.7
% Specific attenuation vs range and total one-way attenuation
clc
clear
close all

%% Given
r_rain = 80000:1:120000;
R_rain = ((-(r_rain/1000 - 100).^2)/8) + 50;
R = [zeros(1, 80000) R_rain]; % mm hr^-1
r_all = 0:1:120000;

lambda = 0.1; % m
theta_e = 1 * pi/180; % rad

%% Calculate Kr for 10 cm wavelength
Kr = 0.000343 * R.^0.97; % dB km^-1

%% Calculate accumulated attenuation - rain
Att_tot = zeros(1, length(r_all));
for i = 2:length(r_all)
   Att_tot(i) = Att_tot(i-1) + Kr(i)*0.001; 
end

%% Calculate gaseous attenuation
Att_gas_2_way = (0.4 + 3.45*exp(-theta_e/1.8)) * ...
    (1 - exp(-(r_all/1000)/(27.8 + 154 * exp(-theta_e/2.2))));
Att_gas = Att_gas_2_way / 2;

%% Plot
% Kr vs r
figure(1)
plot(r_all/1000, Kr)
xlabel('Range [km]')
ylabel('Specific attenuation [dB km^{-1}]')
title('Specific Attenuation vs. Range')
grid on

% Accumulated attenuation
figure(2)
hold on
plot(r_all/1000, Att_tot, 'LineStyle', '--', 'Color', 'k')
plot(r_all/1000, Att_gas, 'LineStyle', ':', 'Color', 'k')
plot(r_all/1000, (Att_tot + Att_gas), 'LineStyle', '-', 'Color', 'k')
grid on
legend('Rain', 'Gas', 'Both')