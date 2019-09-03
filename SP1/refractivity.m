% Brian R. Greene
% refractivity.m
% Plot all data from soundings stored in .m files
clc
clear
close all

%% Load data
% Sounding 1 - KEY Key West Observations at 12Z 03 Sep 2019
KEY_20190903_12Z;

p1 = data(:, 1);
z1 = data(:, 2)/1000;
t1 = data(:, 3);
td1 = data(:, 4);
rh1 = data(:, 5);
r1 = data(:, 6);
wd1 = data(:, 7);
ws1 = data(:, 8);

% Sounding 2 - TUS Tucson Observations at 12Z 03 Sep 2019
TUS_20190903_12Z;

p2 = data(:, 1);
z2 = data(:, 2)/1000;
t2 = data(:, 3);
td2 = data(:, 4);
rh2 = data(:, 5);
r2 = data(:, 6);
wd2 = data(:, 7);
ws2 = data(:, 8);

% Sounding 3 - ILX Lincoln Observations at 12Z 03 Sep 2019
ILX_20190903_12Z;

p3 = data(:, 1);
z3 = data(:, 2)/1000;
t3 = data(:, 3);
td3 = data(:, 4);
rh3 = data(:, 5);
r3 = data(:, 6);
wd3 = data(:, 7);
ws3 = data(:, 8);

%% Calculate mixing ratio and refractivity
% N: refractivity
% Ndry: refractivity w/o water vapor
% Nref: model refractivity (fig 2.7)

% KEY
r1 = r1/1000;
pw1 = p1 .* (r1./(0.622 + r1));
N1 = (77.6./(t1 + 273)) .* (p1 + 4810*(pw1./(t1+273)));
Ndry1 = (77.6./(t1+273)) .* (p1);
Nref1 = 313 * exp(-0.1439*z1);

% TUS
r2 = r2/1000;
pw2 = p2 .* (r2./(0.622 + r2));
N2 = (77.6./(t2 + 273)) .* (p2 + 4810*(pw2./(t2+273)));
Ndry2 = (77.6./(t2+273)) .* (p2);
Nref2 = 313 * exp(-0.1439*z2);

% ILX
r3 = r3/1000;
pw3 = p3 .* (r3./(0.622 + r3));
N3 = (77.6./(t3 + 273)) .* (p3 + 4810*(pw3./(t3+273)));
Ndry3 = (77.6./(t3+273)) .* (p3);
Nref3 = 313 * exp(-0.1439*z3);

%% Plot
f1 = figure(1);
f1.Position = [375,492,1000,600];

% KEY
subplot(1, 3, 1)
hold on
plot(N1, z1, 'DisplayName', 'Refractivity', 'color', 'black')
plot(Ndry1, z1, 'DisplayName', 'Dry Refractivity', 'color', 'red')
plot(Nref1, z1, 'DisplayName', 'Reference Refractivity', 'color', 'blue')
title('Refractivity vs. Altitude for Tropical Location - KEY')
xlabel('Refractivity [-]')
ylabel('Altitude [km]')
ylim([0 18])
legend
grid on

% TUS
subplot(1, 3, 2)
hold on
plot(N2, z2, 'DisplayName', 'Refractivity', 'color', 'black')
plot(Ndry2, z2, 'DisplayName', 'Dry Refractivity', 'color', 'red')
plot(Nref2, z2, 'DisplayName', 'Reference Refractivity', 'color', 'blue')
title('Refractivity vs. Altitude for Desert Location - TUS')
xlabel('Refractivity [-]')
ylabel('Altitude [km]')
legend
grid on

% ILX
subplot(1, 3, 3)
hold on
plot(N3, z3, 'DisplayName', 'Refractivity', 'color', 'black')
plot(Ndry3, z3, 'DisplayName', 'Dry Refractivity', 'color', 'red')
plot(Nref3, z3, 'DisplayName', 'Reference Refractivity', 'color', 'blue')
title('Refractivity vs. Altitude for Inversion - ILX')
xlabel('Refractivity [-]')
ylabel('Altitude [km]')
legend
grid on

% save
print(f1, 'Q6', '-dpng', '-r300')