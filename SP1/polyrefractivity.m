% Brian R. Greene
% polyrefractivity.m
% Plot all data from soundings stored in .m files
clc
clear
close all

%% Load data
% OUN Norman Observations at 12Z 03 Sep 2019
OUN_20190903_12Z;

p = data(:, 1);
z = data(:, 2);
t = data(:, 3);
r = data(:, 6);

%% Calculate refractivity
% N: refractivity
% Ndry: refractivity w/o water vapor
% Nref: model refractivity (fig 2.7)

r = r/1000;
pw = p .* (r./(0.622 + r));
N = (77.6./(t + 273)) .* (p + 4810*(pw./(t+273)));
n = (N * 1e-6) + 1;

%% Fit polynomial
iuse = find(z < 2000);
p = polyfit(z(iuse), n(iuse), 1);
x = 0:1:2000;
fit = p(1)*x + p(2);

%% Find ae
% slope = dN/dz = p(1)
a = 6374e3;
dNdz = p(1);
ae = 1/((1/a) + dNdz);
disp('ae = ');
disp(ae);

%% Plot
f1 = figure(1);
f1.Position = [375,492,600,600];

hold on
scatter(n(iuse), z(iuse))
plot(fit, x)
title('Refractivity n vs. Altitude with polyfit')
xlabel('Refractivity n [-]')
ylabel('Altitude [m]')
grid on

% save
print(f1, 'Q7', '-dpng', '-r300')