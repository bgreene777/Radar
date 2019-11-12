% Brian R. Greene
% Use X_h and X_v to estimate differential reflectivity ZDR

clc
clear
close all

%% Load data
load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200821_e02.60.mat');

%% Get actual r values
r_all = linspace(0, 1009, 1010) * delr / 1000;

%% Calc power
X_h = double(X_h);
X_v = double(X_v);
el_rad = el/180*pi;

[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

% all pulses
R_h=mean(X_h(:,:,1:num_pulses).*conj(X_h(:,:,1:num_pulses)),3);
R_v=mean(X_v(:,:,1:num_pulses).*conj(X_v(:,:,1:num_pulses)),3);

% calc ZDR
Zdr = 10*log10(R_h ./ R_v);

%% Plot

f1=figure(1);
set(gcf,'render','painters', 'Position', [10 10 900 900]);
% all pulses
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,Zdr);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[Zdr;Zdr(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('carbmap'));
caxis([-3 8]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{DR} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c = colorbar;
c.Label.String = 'Differential Reflectivity (dB)';
% save
print(f1, 'Q4_Z_DR', '-dpng', '-r300')