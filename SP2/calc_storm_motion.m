% Brian R. Greene
% calculate storm motion based off difference in position of centroid of
% top 5% of horizontal power return locations

clc
clear
close all

%% Load data
dat1 = load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200821_e02.60.mat');
dat2 = load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200841_e02.60.mat');

%% Get actual r values
r_all = linspace(0, 1009, 1010) * dat1.delr / 1000;
%% Calculate powers
% frame 1
dat1.X_h = double(dat1.X_h);
el_rad = dat1.el/180*pi;

[r,az_rad] = meshgrid(((0:dat1.num_gates-1)*dat1.delr+dat1.r_min)/1e3,dat1.az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad);

dat1.R_h=mean(dat1.X_h(:,:,1:dat1.num_pulses).*conj(dat1.X_h(:,:,1:dat1.num_pulses)),3);

dat1.pow_h = 10*log10(dat1.R_h);

% frame 2
dat2.X_h = double(dat2.X_h);

dat2.R_h=mean(dat2.X_h(:,:,1:dat2.num_pulses).*conj(dat2.X_h(:,:,1:dat2.num_pulses)),3);

dat2.pow_h = 10*log10(dat2.R_h);

%% Calculate 95th percentile of powers in each frame
% frame 1
p95_1 = prctile(dat1.pow_h, 95, [1 2]);

[iaz95_1, ir95_1] = find(dat1.pow_h >= p95_1);

r95_1 = r_all(ir95_1);
az95_1 = dat1.az_set(iaz95_1);
x95_1 = r95_1 * cos(el_rad) .* sind(az95_1);
y95_1 = r95_1 * cos(el_rad) .* cosd(az95_1);

x95_1_mean = median(x95_1);
y95_1_mean = median(y95_1);

% frame 2
p95_2 = prctile(dat2.pow_h, 95, [1 2]);

[iaz95_2, ir95_2] = find(dat2.pow_h >= p95_2);

r95_2 = r_all(ir95_2);
az95_2 = dat2.az_set(iaz95_2);
x95_2 = r95_2 * cos(el_rad) .* sind(az95_2);
y95_2 = r95_2 * cos(el_rad) .* cosd(az95_2);

x95_2_mean = median(x95_2);
y95_2_mean = median(y95_2);

%% Calc storm motion
dx = x95_2_mean - x95_1_mean;
dy = y95_2_mean - y95_1_mean;
dd = sqrt(dx^2 + dy^2);
% speed in km/hr
spd = 3600 * dd / 20;
direction = atan2d(dy,dx);

%% Plot
% frame 1
f1=figure(1); 
set(gcf,'render','painters');
naz_max = 360;
if length(dat1.az_set)<naz_max
	pcolor(x,y,dat1.pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dat1.pow_h;dat1.pow_h(1,:)]);
end
hold on
plot(x95_1_mean, y95_1_mean, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(dat1.scan_time),' El=',num2str(dat1.el,'%5.2f'),' degrees']);
colorbar;
% save
print(f1, 'Q8_frame1', '-dpng', '-r300')

% frame 2
f2=figure(2); 
set(gcf,'render','painters');
if length(dat2.az_set)<naz_max
	pcolor(x,y,dat2.pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dat2.pow_h;dat2.pow_h(1,:)]);
end
hold on
plot(x95_2_mean, y95_2_mean, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(dat2.scan_time),' El=',num2str(dat2.el,'%5.2f'),' degrees']);
colorbar;
% save
print(f2, 'Q8_frame2', '-dpng', '-r300')