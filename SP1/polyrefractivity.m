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

%% Calculate gradient
iuse = find(z < 2000);
dNdz = nan(length(iuse)-1, 1);
znew = nan(length(iuse)-1, 1);

for i = 2:length(iuse)
    dN = N(i) - N(i-1);
    dz = z(i) - z(i-1);
    dNdz(i-1) = dN / dz;
    znew(i-1) = (z(i) + z(i-1)) / 2;
end

%% Fit polynomial
p = polyfit(znew, dNdz, 1);
x = 0:1:2000;
fit = p(1)*x + p(2);

%% Plot
f1 = figure(1);
f1.Position = [375,492,600,600];

hold on
scatter(dNdz, znew)
plot(fit, x)