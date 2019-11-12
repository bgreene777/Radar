% Brian R. Greene
% Use X_h and X_v to estimate cross-correlation coefficient rho hv

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

% loop to calculate rhohv
rho_hv = zeros(num_az, num_gates);
for i = 1:num_az
    for j = 1:num_gates
        % all pulses
        [Rhv, ~] = xcorr(squeeze(X_h(i,j,1:num_pulses)),...
            squeeze(X_v(i,j,1:num_pulses)), 0, 'biased');
        rho_hv(i,j) = abs(Rhv) / sqrt(R_h(i,j) * R_v(i,j));
    end
end
disp('done')

%% Plot

f1=figure(1);
set(gcf,'render','painters', 'Position', [10 10 900 900]);
% all pulses
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,rho_hv);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[rho_hv;rho_hv(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('carbmap'));
caxis([0.7 1]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} \rho_{HV} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c = colorbar;
c.Label.String = 'Correlation Coefficient';
% save
print(f1, 'Q5_Rho_HV', '-dpng', '-r300')