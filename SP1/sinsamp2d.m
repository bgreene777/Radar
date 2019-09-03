% Brian R. Greene
% 2dsinsamp.m
%
% This program will sample a 2d sinusoidal signal
% and plot the result using mesh, contour, contourf
%
clc
clear
close all

% samp = 10 Hz
samptime = 1/10;
x = [-20:20] * samptime;
y = [-20:20] * samptime;

[xx, yy] = meshgrid(x, y);

sig = sin(2*xx + 3*yy);

%% Plot
f1 = figure(1);
f1.Position = [797,412,1260,874];

% mesh
subplot(2,2,1)
mesh(xx, yy, sig)
ax1 = gca;
ax1.FontSize = 14;
xlabel('x', 'fontsize', 20)
ylabel('y', 'fontsize', 20)
zlabel('signal', 'fontsize', 20)
title('(a) sin(2x+3y) using mesh', 'fontsize', 20)

% contour
subplot(2,2,2)
contour(xx, yy, sig)
ax1 = gca;
ax1.FontSize = 14;
xlabel('x', 'fontsize', 20)
ylabel('y', 'fontsize', 20)
title('(b) sin(2x+3y) using contour', 'fontsize', 20)

% contourf
subplot(2,2,3)
contourf(xx, yy, sig)
ax1 = gca;
ax1.FontSize = 14;
xlabel('x', 'fontsize', 20)
ylabel('y', 'fontsize', 20)
title('(c) sin(2x+3y) using contourf', 'fontsize', 20)

% waterfall
subplot(2,2,4)
waterfall(xx, yy, sig)
ax1 = gca;
ax1.FontSize = 14;
xlabel('x', 'fontsize', 20)
ylabel('y', 'fontsize', 20)
zlabel('signal', 'fontsize', 20)
title('(d) sin(2x+3y) using waterfall', 'fontsize', 20)

% save
print(f1, 'Q3', '-dpng', '-r200')