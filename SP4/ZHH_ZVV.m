% Brian R. Greene
% Use X_h and X_v to estimate horizontal and vertical reflectivity factors

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
pow_h = 10*log10(R_h);
pow_v = 10*log10(R_v);
% 5 pulses
R_h_5 = mean(X_h(:,:,1:5).*conj(X_h(:,:,1:5)),3);
R_v_5 = mean(X_v(:,:,1:5).*conj(X_v(:,:,1:5)),3);
pow_h_5 = 10*log10(R_h_5);
pow_v_5 = 10*log10(R_v_5);
% 1 pulse
R_h_1 = X_h(:,:,1).*conj(X_h(:,:,1));
R_v_1 = X_v(:,:,1).*conj(X_v(:,:,1));
pow_h_1 = 10*log10(R_h_1);
pow_v_1 = 10*log10(R_v_1);
% correct for range
r_cal_dB = 20*log10(1e3*r);
% calib to same range as KTLX
pow_cal_dB = -80;

% calc dBZ
dBZ_h = pow_h + r_cal_dB + pow_cal_dB;
dBZ_v = pow_v + r_cal_dB + pow_cal_dB;

dBZ_h_5 = pow_h_5 + r_cal_dB + pow_cal_dB;
dBZ_v_5 = pow_v_5 + r_cal_dB + pow_cal_dB;

dBZ_h_1 = pow_h_1 + r_cal_dB + pow_cal_dB;
dBZ_v_1 = pow_v_1 + r_cal_dB + pow_cal_dB;

%% Z_HH PPI - 50, 5, and 1 pulse

f1=figure(1);
set(gcf,'render','painters', 'Position', [10 10 1800 600]);
% all pulses
subplot(1, 3, 1)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,dBZ_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_h;dBZ_h(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{HH} all ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% 5 pulses
subplot(1, 3, 2)
if length(az_set)<naz_max
	pcolor(x,y,dBZ_h_5);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_h_5;dBZ_h_5(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{HH} 5 pulses ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% 1 pulse
subplot(1, 3, 3)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,dBZ_h_1);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_h_1;dBZ_h_1(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{HH} 1 pulse ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f1, 'Q1_Z_HH', '-dpng', '-r300')

%% Z_VV PPI

f2=figure(2); 
set(gcf,'render','painters', 'Position', [10 10 1800 600]);
% all pulses
subplot(1, 3, 1)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,dBZ_v);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_v;dBZ_v(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{VV} all ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% 5 pulses
subplot(1, 3, 2)
if length(az_set)<naz_max
	pcolor(x,y,dBZ_v_5);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_v_5;dBZ_v_5(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{VV} 5 pulses ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% 1 pulse
subplot(1, 3, 3)
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,dBZ_v_1);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ_v_1;dBZ_v_1(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
caxis([10 80]);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Z_{VV} 1 pulse ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f2, 'Q1_Z_VV', '-dpng', '-r300')