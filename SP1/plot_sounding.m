% Brian R. Greene
% plot_sounding.m
% Plot all data from soundings stored in .m files
clc
clear
close all

%% Load data
% Sounding 1 - OUN Norman Observations at 12Z 03 Sep 2019
OUN_20190903_12Z;

p1 = data(:, 1);
z1 = data(:, 2);
t1 = data(:, 3);
td1 = data(:, 4);
rh1 = data(:, 5);
r1 = data(:, 6);
wd1 = data(:, 7);
ws1 = data(:, 8);
th1 = data(:, 9);
the1 = data(:, 10);
thv1 = data(:, 11);

% Sounding 2 - KEY Key West Observations at 12Z 03 Sep 2019
KEY_20190903_12Z;

p2 = data(:, 1);
z2 = data(:, 2);
t2 = data(:, 3);
td2 = data(:, 4);
rh2 = data(:, 5);
r2 = data(:, 6);
wd2 = data(:, 7);
ws2 = data(:, 8);
th2 = data(:, 9);
the2 = data(:, 10);
thv2 = data(:, 11);

%% Plot
% Sounding 1
f1 = figure(1);
f1.Position = [375,492,1200,1400];

% Temperature, Dewpoint
subplot(4, 2, 1)
hold on
plot(t1, z1, 'color', 'r', 'DisplayName', 'Temperature')
plot(td1, z1, 'color', 'g', 'DisplayName', 'Dewpoint')
legend
title('Temperature and Dewpoint vs. Altitude')
xlabel('Temperature and Dewpoint Temperature [C]')
ylabel('Altitude [m]')
grid on

% RH
subplot(4, 2, 2)
plot(rh1, z1, 'color', 'b')
title('Relative Humidity vs. Altitude')
xlabel('Relative Humidity [percent]')
ylabel('Altitude [m]')
grid on

% mixing ratio
subplot(4, 2, 3)
plot(r1, z1, 'color', 'b')
title('Mixing Ratio vs. Altitude')
xlabel('Mixing Ratio [g/kg]')
ylabel('Altitude [m]')
grid on

% wind direction
subplot(4, 2, 4)
plot(wd1, z1, 'color', 'b')
title('Wind Direction vs. Altitude')
xlabel('Wind Direction [deg]')
ylabel('Altitude [m]')
grid on

% wind speed
subplot(4, 2, 5)
plot(ws1, z1, 'color', 'b')
title('Wind Speed vs. Altitude')
xlabel('Wint Speed [kts]')
ylabel('Altitude [m]')
grid on

% Theta, Theta-e, Theta-v
subplot(4, 2, 6)
hold on
plot(th1, z1, 'color', 'k', 'DisplayName', 'Theta')
plot(the1, z1, 'color', 'b', 'DisplayName', 'Theta-e')
plot(thv1, z1, 'color', 'r', 'DisplayName', 'Theta-v')
legend
title('Potential Temperature vs. Altitude')
xlabel('Potential Temperature [K]')
ylabel('Altitude [m]')
grid on

% pressure
subplot(4, 2, 7)
plot(p1, z1, 'color', 'b')
title('Pressure vs. Altitude')
xlabel('Pressure [hPa]')
ylabel('Altitude [m]')
grid on

% save
print(f1, 'Q5_OUN', '-dpng', '-r300')

% Sounding 2
f2 = figure(2);
f2.Position = [375,492,1200,1400];

% Temperature, Dewpoint
subplot(4, 2, 1)
hold on
plot(t2, z2, 'color', 'r', 'DisplayName', 'Temperature')
plot(td2, z2, 'color', 'g', 'DisplayName', 'Dewpoint')
legend
title('Temperature and Dewpoint vs. Altitude')
xlabel('Temperature and Dewpoint Temperature [C]')
ylabel('Altitude [m]')
grid on

% RH
subplot(4, 2, 2)
plot(rh2, z2, 'color', 'b')
title('Relative Humidity vs. Altitude')
xlabel('Relative Humidity [percent]')
ylabel('Altitude [m]')
grid on

% mixing ratio
subplot(4, 2, 3)
plot(r2, z2, 'color', 'b')
title('Mixing Ratio vs. Altitude')
xlabel('Mixing Ratio [g/kg]')
ylabel('Altitude [m]')
grid on

% wind direction
subplot(4, 2, 4)
plot(wd2, z2, 'color', 'b')
title('Wind Direction vs. Altitude')
xlabel('Wind Direction [deg]')
ylabel('Altitude [m]')
grid on

% wind speed
subplot(4, 2, 5)
plot(ws2, z2, 'color', 'b')
title('Wind Speed vs. Altitude')
xlabel('Wint Speed [kts]')
ylabel('Altitude [m]')
grid on

% Theta, Theta-e, Theta-v
subplot(4, 2, 6)
hold on
plot(th2, z2, 'color', 'k', 'DisplayName', 'Theta')
plot(the2, z2, 'color', 'b', 'DisplayName', 'Theta-e')
plot(thv2, z2, 'color', 'r', 'DisplayName', 'Theta-v')
legend
title('Potential Temperature vs. Altitude')
xlabel('Potential Temperature [K]')
ylabel('Altitude [m]')
grid on

% pressure
subplot(4, 2, 7)
plot(p2, z2, 'color', 'b')
title('Pressure vs. Altitude')
xlabel('Pressure [hPa]')
ylabel('Altitude [m]')
grid on

% save
print(f2, 'Q5_KEY', '-dpng', '-r300')