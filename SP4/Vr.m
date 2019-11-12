% Brian R. Greene
% pulse-pair estimator for radial velocity

clc
clear
close all

%% Load data
load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200821_e02.60.mat');

%% Get actual r values
r_all = linspace(0, 1009, 1010) * delr / 1000;

%% Calc vr
X_h = double(X_h);
X_v = double(X_v);
el_rad = el/180*pi;

[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

% use ppp
vr_h = zeros(num_az, num_gates);
vr_v = zeros(num_az, num_gates);
vr_h_5 = zeros(num_az, num_gates);
vr_v_5 = zeros(num_az, num_gates);
coeff = -lambda/(4*pi*pri);
for i = 1:num_az
    for j = 1:num_gates
        % all pulses
        [Ri_h, ~] = xcorr(squeeze(X_h(i,j,1:num_pulses)), 1, 'biased');
        vr_h(i,j) = coeff * angle(Ri_h(1)); % take lag = 1
        [Ri_v, ~] = xcorr(squeeze(X_v(i,j,1:num_pulses)), 1, 'biased');
        vr_v(i,j) = coeff * angle(Ri_v(1)); % take lag = 1
        % 5 pulses
        [Ri_h_5, ~] = xcorr(squeeze(X_h(i,j,1:5)), 1, 'biased');
        vr_h_5(i,j) = coeff * angle(Ri_h_5(1)); % take lag = 1
        [Ri_v_5, ~] = xcorr(squeeze(X_v(i,j,1:5)), 1, 'biased');
        vr_v_5(i,j) = coeff * angle(Ri_v_5(1)); % take lag = 1
    end
end
disp('done')

%% Vr_h

f1=figure(1);
set(gcf,'render','painters', 'Position', [10 10 1200 600]);
% all pulses
subplot(1, 2, 1)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,vr_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[vr_h;vr_h(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('brmap'));
caxis([-20 20])
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} V_{r, horizontal}^{all} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c1 = colorbar;
c1.Label.String = 'Radial Velocity (m/s)';
% 5 pulses
subplot(1, 2, 2)
if length(az_set)<naz_max
	pcolor(x,y,vr_h_5);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[vr_h_5;vr_h_5(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('brmap'));
caxis([-20 20])
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} V_{r, horizontal}^{5 pulses} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c2 = colorbar;
c2.Label.String = 'Radial Velocity (m/s)';
% save
print(f1, 'Q3_vr_h', '-dpng', '-r300')

%% Vr_v

f2=figure(2);
set(gcf,'render','painters', 'Position', [10 10 1200 600]);
% all pulses
subplot(1, 2, 1)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,vr_v);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[vr_v;vr_v(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('brmap'));
caxis([-20 20])
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} V_{r, vertical}^{all} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c1 = colorbar;
c1.Label.String = 'Radial Velocity (m/s)';
% 5 pulses
subplot(1, 2, 2)
if length(az_set)<naz_max
	pcolor(x,y,vr_v_5);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[vr_v_5;vr_v_5(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(boonlib('brmap'));
caxis([-20 20])
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} V_{r, vertical}^{5 pulses} ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
c2 = colorbar;
c2.Label.String = 'Radial Velocity (m/s)';
% save
print(f2, 'Q3_vr_v', '-dpng', '-r300')