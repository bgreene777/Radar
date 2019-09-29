% Brian R. Greene
% Plot I(t) and Q(t) together for several aximuth angles and range gates

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

R_h=mean(X_h(:,:,1:num_pulses).*conj(X_h(:,:,1:num_pulses)),3);
R_v=mean(X_v(:,:,1:num_pulses).*conj(X_v(:,:,1:num_pulses)),3);

pow_h = 10*log10(R_h);
pow_v = 10*log10(R_v);

%% Inside storm: range gate 600, azimuth angle 2
r_1 = r_all(600);
az_1 = az_set(2);
x_1 = r_1 * cos(el_rad).*sind(az_1);
y_1 = r_1 * cos(el_rad).*cosd(az_1);

% Position in PPI
f1=figure(1); 
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_h;pow_h(1,:)]);
end
hold on
plot(x_1, y_1, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f1, 'Q6_storm_PPI', '-dpng', '-r300')

% X_h
f2=figure(2);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_h(2,600,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_h(2,600,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_h at az = ', num2str(az_1, '%3.0f'), '^\circ, r = ', num2str(r_1, '%6.3f'), ' km'])
% save
print(f2, 'Q6_storm_h', '-dpng', '-r300')

% X_v
f3=figure(3);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_v(2,600,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_v(2,600,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_v at az = ', num2str(az_1, '%3.0f'), '^\circ, r = ', num2str(r_1, '%6.3f'), ' km'])
% save
print(f3, 'Q6_storm_v', '-dpng', '-r300')

%% Close to radar: range gate 15, azimuth angle 270
r_2 = r_all(15);
az_2 = az_set(270);
x_2 = r_2 * cos(el_rad).*sind(az_2);
y_2 = r_2 * cos(el_rad).*cosd(az_2);

% Position in PPI
f4=figure(4); 
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_h;pow_h(1,:)]);
end
hold on
plot(x_2, y_2, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f4, 'Q6_clutter_PPI', '-dpng', '-r300')

% X_h
f5=figure(5);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_h(270,15,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_h(270,15,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_h at az = ', num2str(az_2, '%3.0f'), '^\circ, r = ', num2str(r_2, '%6.3f'), ' km'])
% save
print(f5, 'Q6_clutter_h', '-dpng', '-r300')

% X_v
f6=figure(6);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_v(270,15,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_v(270,15,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_v at az = ', num2str(az_2, '%3.0f'), '^\circ, r = ', num2str(r_2, '%6.3f'), ' km'])
% save
print(f6, 'Q6_clutter_v', '-dpng', '-r300')

%% Far from radar: range gate 900, azimuth angle 180
r_3 = r_all(900);
az_3 = az_set(180);
x_3 = r_3 * cos(el_rad).*sind(az_3);
y_3 = r_3 * cos(el_rad).*cosd(az_3);

% Position in PPI
f7=figure(7); 
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_h;pow_h(1,:)]);
end
hold on
plot(x_3, y_3, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f7, 'Q6_clear_PPI', '-dpng', '-r300')

% X_h
f8=figure(8);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_h(180,900,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_h(180,900,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_h at az = ', num2str(az_3, '%3.0f'), '^\circ, r = ', num2str(r_3, '%6.3f'), ' km'])
% save
print(f8, 'Q6_clear_h', '-dpng', '-r300')

% X_v
f9=figure(9);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_v(180,900,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_v(180,900,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_v at az = ', num2str(az_3, '%3.0f'), '^\circ, r = ', num2str(r_3, '%6.3f'), ' km'])
% save
print(f9, 'Q6_clear_v', '-dpng', '-r300')
