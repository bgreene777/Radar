% Brian R. Greene
% Plot I(t) and Q(t) together for several aximuth angles and range gates
% Calculate autocorrelation function for each timeseries

clc
clear
close all

%% Load data
load('/Users/briangreene/Documents/MATLAB/Radar/SP2/Data/iq_PX-1000_20130520_200821_e02.60.mat');

%% Get actual r values
r_all = linspace(0, 1009, 1010) * delr / 1000;

%% Calc power
X_h = double(X_h);
el_rad = el/180*pi;

[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

R_h=mean(X_h(:,:,1:num_pulses).*conj(X_h(:,:,1:num_pulses)),3);

pow_h = 10*log10(R_h);

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
print(f1, 'Q1_storm_PPI', '-dpng', '-r300')

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
print(f2, 'Q1_storm_h', '-dpng', '-r300')

% ACF
[R_1, lags_1] = xcorr(squeeze(X_h(2,600,:)));
lags_1 = lags_1 * pri;
f3=figure(3);
plot(lags_1, R_1)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Autocorrelation Function')
title('\bf \fontsize{12} Autocorrelation for storm data')
% save
print(f3, 'Q1_storm_acf', '-dpng', '-r300')

% normalized ACF
R_1_norm =  R_1/max(R_1);
id_1 = find(R_1_norm(50:end) <= 0.7);
t_dcor_1 = lags_1(50+id_1(1)-2);
R_dcor_1 = R_1_norm(50+id_1(1)-2);
f31 = figure(31);
hold on
plot(lags_1, R_1_norm)
plot(lags_1, 0.7*ones(1,length(lags_1)))
plot(t_dcor_1, R_dcor_1, 'r*', 'MarkerSize', 16)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Normalized Autocorrelation Function')
title(['\bf \fontsize{12} Normalized Autocorrelation for storm data, Decorrelation time=',num2str(t_dcor_1, '%5.5f'), ' s'])
% save
print(f31, 'Q2_storm_acf_norm', '-dpng', '-r300')

%% Clutter: range gate 15, azimuth angle 270
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
print(f4, 'Q1_clutter_PPI', '-dpng', '-r300')

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
print(f5, 'Q1_clutter_h', '-dpng', '-r300')

% ACF
[R_2, lags_2] = xcorr(squeeze(X_h(270,15,:)));
lags_2 = lags_2 * pri;
f6=figure(6);
plot(lags_2, R_2)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Autocorrelation Function')
title('\bf \fontsize{12} Autocorrelation for clutter')
% save
print(f6, 'Q1_clutter_acf', '-dpng', '-r300')

% normalized ACF
R_2_norm =  R_2/max(R_2);
id_2 = find(R_2_norm(50:end) <= 0.7);
t_dcor_2 = lags_2(50+id_2(1)-2);
R_dcor_2 = R_2_norm(50+id_2(1)-2);
f61 = figure(61);
hold on
plot(lags_2, R_2_norm)
plot(lags_2, 0.7*ones(1,length(lags_1)))
plot(t_dcor_2, R_dcor_2, 'r*', 'MarkerSize', 16)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Normalized Autocorrelation Function')
title(['\bf \fontsize{12} Normalized Autocorrelation for ground clutter, Decorrelation time=',num2str(t_dcor_2, '%5.5f'), ' s'])
% save
print(f61, 'Q2_clutter_acf_norm', '-dpng', '-r300')
%% Clear air: range gate 900, azimuth angle 180
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
print(f7, 'Q1_clear_PPI', '-dpng', '-r300')

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
print(f8, 'Q1_clear_h', '-dpng', '-r300')

% ACF
[R_3, lags_3] = xcorr(squeeze(X_h(180,900,:)));
lags_3 = lags_3 * pri;
f9=figure(9);
plot(lags_3, R_3)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Autocorrelation Function')
title('\bf \fontsize{12} Autocorrelation for clear air')
% save
print(f9, 'Q1_clear_acf', '-dpng', '-r300')

% normalized ACF
R_3_norm =  R_3/max(R_3);
id_3 = find(R_3_norm(50:end) <= 0.7);
t_dcor_3 = lags_1(50+id_3(1)-2);
R_dcor_3 = R_1_norm(50+id_3(1)-2);
f91 = figure(91);
hold on
plot(lags_3, R_3_norm)
plot(lags_3, 0.7*ones(1,length(lags_3)))
plot(t_dcor_3, R_dcor_3, 'r*', 'MarkerSize', 16)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Normalized Autocorrelation Function')
title(['\bf \fontsize{12} Normalized Autocorrelation for clear air, Decorrelation time=',num2str(t_dcor_3, '%5.5f'), ' s'])
% save
print(f91, 'Q2_clear_acf_norm', '-dpng', '-r300')
%% Tornado: range gate 437, azimuth angle 53
r_4 = r_all(437);
az_4 = az_set(53);
x_4 = r_4 * cos(el_rad).*sind(az_4);
y_4 = r_4 * cos(el_rad).*cosd(az_4);

f10=figure(10); 
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_h;pow_h(1,:)]);
end
hold on
plot(x_4, y_4, 'r*', 'MarkerSize', 16)
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f10, 'Q1_tor_PPI', '-dpng', '-r300')

% X_h
f11=figure(11);
hold on
plot([1:num_pulses]*pri, real(squeeze(X_h(53,437,:))), 'DisplayName', 'I(t)')
plot([1:num_pulses]*pri, imag(squeeze(X_h(53,437,:))), 'DisplayName', 'Q(t)')
grid on
legend('show');
xlabel('\bf \fontsize{11} Time (s)')
ylabel('\bf \fontsize{11} I and Q')
title(['\bf \fontsize{12} X_h at az = ', num2str(az_4, '%3.0f'), '^\circ, r = ', num2str(r_4, '%6.3f'), ' km'])
% save
print(f11, 'Q1_tor_h', '-dpng', '-r300')

% ACF
[R_4, lags_4] = xcorr(squeeze(X_h(53,437,:)));
lags_4 = lags_4 * pri;
f12=figure(12);
plot(lags_4, R_4)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Autocorrelation Function')
title('\bf \fontsize{12} Autocorrelation for tornado data')
% save
print(f12, 'Q1_tor_acf', '-dpng', '-r300')

% normalized ACF
R_4_norm =  R_4/max(R_4);
id_4 = find(R_4_norm(50:end) <= 0.7);
t_dcor_4 = lags_4(50+id_4(1)-2);
R_dcor_4 = R_4_norm(50+id_4(1)-2);
f121 = figure(121);
hold on
plot(lags_4, R_4_norm)
plot(lags_4, 0.7*ones(1,length(lags_4)))
plot(t_dcor_4, R_dcor_4, 'r*', 'MarkerSize', 16)
grid on
xlabel('\bf \fontsize{11} Lag (s)')
ylabel('\bf \fontsize{11} Normalized Autocorrelation Function')
title(['\bf \fontsize{12} Normalized Autocorrelation for tornado, Decorrelation time=',num2str(t_dcor_4, '%5.5f'), ' s'])
% save
print(f121, 'Q2_tor_acf_norm', '-dpng', '-r300')

%% Scatter plots: storm data - range gate 600, azimuth angle 2
f30 = figure(30);
hold on
plot(real(squeeze(X_h(2, 600, :))), imag(squeeze(X_h(2, 600, :))), 'b.')
plot(real(squeeze(X_h(2, 599, :))), imag(squeeze(X_h(2, 599, :))), 'b.')
plot(real(squeeze(X_h(2, 601, :))), imag(squeeze(X_h(2, 601, :))), 'b.')
plot(real(squeeze(X_h(1, 600, :))), imag(squeeze(X_h(1, 600, :))), 'b.')
plot(real(squeeze(X_h(3, 600, :))), imag(squeeze(X_h(3, 600, :))), 'b.')
axis([-1000 1000 -1000 1000])
axis square
grid on
xlabel('\bf \fontsize{11} I(t)')
ylabel('\bf \fontsize{11} Q(t)')
title('\bf \fontsize{12} I/Q scatter for data inside storm')
% save
print(f30, 'Q3_iq_storm', '-dpng', '-r300')

%% Scatter plots: clear air
v_2 = cat(1, squeeze(X_h(180, 900, :)), squeeze(X_h(180, 899, :)), squeeze(X_h(180, 901, :)), squeeze(X_h(181, 900, :)), squeeze(X_h(179, 900, :)));
I_2 = real(v_2);
Q_2 = imag(v_2);

f35 = figure(35);
plot(I_2, Q_2, 'b.')

axis([-10 10 -10 10])
axis square
grid on
xlabel('\bf \fontsize{11} I(t)')
ylabel('\bf \fontsize{11} Q(t)')
title('\bf \fontsize{12} I/Q scatter for data clear air')
% save
print(f35, 'Q3_iq_clear', '-dpng', '-r300')

%% Histogram: storm data - range gate 600, azimuth angle 2
I = cat(1, real(squeeze(X_h(2, 600, :))), real(squeeze(X_h(2, 599, :))), real(squeeze(X_h(2, 601, :))), real(squeeze(X_h(1, 600, :))), real(squeeze(X_h(3, 600, :))));
Q = cat(1, imag(squeeze(X_h(2, 600, :))), imag(squeeze(X_h(2, 599, :))), imag(squeeze(X_h(2, 601, :))), imag(squeeze(X_h(1, 600, :))), imag(squeeze(X_h(3, 600, :))));

f31 = figure(31);
histogram(I, -1000:100:1000)
xlabel('\bf \fontsize{11} I(t)')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} I(t) distribution inside storm')
% save
print(f31, 'Q4_i_hist', '-dpng', '-r300')
f32 = figure(32);
histogram(Q, -1000:100:1000)
xlabel('\bf \fontsize{11} Q(t)')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Q(t) distribution inside storm')
% save
print(f32, 'Q4_q_hist', '-dpng', '-r300')

%% Histogram: clear air, range 900 az 180
f33 = figure(33);
histogram(I_2)
xlabel('\bf \fontsize{11} I(t)')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} I(t) distribution clear air')
% save
print(f33, 'Q4_i_hist_clear', '-dpng', '-r300')

f34 = figure(34);
histogram(Q)
xlabel('\bf \fontsize{11} Q(t)')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Q(t) distribution clear air')
% save
print(f34, 'Q4_q_hist_clear', '-dpng', '-r300')

%% Amplitude, phase, power
v = cat(1, squeeze(X_h(2, 600, :)), squeeze(X_h(2, 599, :)), squeeze(X_h(2, 601, :)), squeeze(X_h(1, 600, :)), squeeze(X_h(3, 600, :)));

% amplitude
v_amp = sqrt(v.*conj(v));
f40 = figure(40);
histogram(v_amp, 0:50:1000)
xlabel('\bf \fontsize{11} Amplitude of echo voltage')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage amplitude distribution inside storm')
% save
print(f40, 'Q5_amp_hist', '-dpng', '-r300')

% phase
v_phase = atan2d(Q, I);
f41 = figure(41);
histogram(v_phase, -200:20:200)
xlabel('\bf \fontsize{11} Echo voltage phase')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage phase distribution inside storm')
% save
print(f41, 'Q5_phase_hist', '-dpng', '-r300')

% power
v_pow = v.*conj(v);
f42 = figure(42);
histogram(v_pow)
xlabel('\bf \fontsize{11} Echo voltage power')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage power distribution inside storm')
% save
print(f42, 'Q5_pow_hist', '-dpng', '-r300')

%% Amplitude, phase, power clear air

% amplitude
v_amp_2 = sqrt(v_2.*conj(v_2));
f43 = figure(43);
histogram(v_amp_2)
xlabel('\bf \fontsize{11} Amplitude of echo voltage')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage amplitude distribution clear air')
% save
print(f43, 'Q5_amp_hist_clear', '-dpng', '-r300')

% phase
v_phase_2 = atan2d(imag(v_2), real(v_2));
f44 = figure(44);
histogram(v_phase_2, -200:20:200)
xlabel('\bf \fontsize{11} Echo voltage phase')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage phase distribution clear air')
% save
print(f44, 'Q5_phase_hist_clear', '-dpng', '-r300')

% power
v_pow_2 = v_2.*conj(v_2);
f45 = figure(45);
histogram(v_pow_2)
xlabel('\bf \fontsize{11} Echo voltage power')
ylabel('\bf \fontsize{11} Count')
title('\bf \fontsize{12} Echo voltage power distribution clear air')
% save
print(f45, 'Q5_pow_hist_clear', '-dpng', '-r300')

%% FFT - storm
Fs = 1./pri;
va = lambda * Fs / 4;
L = 2^nextpow2(length(R_1_norm));
S_1 = fft(R_1_norm, L);
S_1_shift = abs(fftshift(S_1));
n = length(S_1);
f0 = (0:n-1) * (2 * va/n);
f = (-n/2:n/2-1) * (2 * va/n);
f60 = figure(60);
subplot(2, 1, 1)
plot(f0, 10*log10(S_1))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Power Spectral Density In Storm')
subplot(2, 1, 2)
plot(f, 10*log10(S_1_shift))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Shifted Power Spectral Density In Storm')
% subplot(3,1,3)
% plot(10*log10(pow))
% save
print(f60, 'Q6_psd_fft_storm', '-dpng', '-r300')

%% FFT - clutter
S_2 = fft(R_2_norm, L);
S_2_shift = fftshift(S_2);
f61 = figure(61);
subplot(2,1,1)
plot(f0, 10*log10(S_2))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Power Spectral Density In Clutter')
subplot(2,1,2)
plot(f, 10*log10(S_2_shift))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Shifted Power Spectral Density In Clutter')
% save
print(f61, 'Q6_psd_fft_clutter', '-dpng', '-r300')

%% FFT - clear air
S_3 = fft(R_3_norm, L);
S_3_shift = fftshift(S_3);
f62 = figure(62);
subplot(2,1,1)
plot(f0, 10*log10(S_3))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Power Spectral Density In Clear Air')
subplot(2,1,2)
plot(f, 10*log10(S_3_shift))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Shifted Power Spectral Density In Clear Air')
% save
print(f62, 'Q6_psd_fft_storm_clear', '-dpng', '-r300')

%% FFT - tornado
S_4 = fft(R_4_norm, L);
S_4_shift = fftshift(S_4);
f63 = figure(63);
subplot(2,1,1)
plot(f0, 10*log10(S_4))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Power Spectral Density In Tornado')
subplot(2,1,2)
plot(f, 10*log10(S_4_shift))
xlabel('\bf \fontsize{11} Doppler Velocity (m/s)')
ylabel('\bf \fontsize{11} Power Spectral Density (dB)')
title('\bf \fontsize{12} Shifted Power Spectral Density In Tornado')
% save
print(f63, 'Q6_psd_fft_tor', '-dpng', '-r300')

%% All ACF; integrate PSD
R_all = zeros(num_az, num_gates);
Pow_all = zeros(num_az, num_gates);
f_all = (-n/2:n/2-1)*Fs/n;
for i = 1:num_az
    for j = 1:num_gates
        [Ri, ~] = xcorr(squeeze(X_h(i,j,1:num_pulses)));
        R_all(i,j) = real(Ri(50)); % take middle element (0th index)
        Si = abs(fftshift(fft(Ri, L)));
        Pow_all(i,j) = trapz(f_all, Si);
    end
end
R_all_dB = 10*log10(R_all);
Pow_all_dB = 10*log10(Pow_all);
disp('done')

%% PPI - xcorr
f70=figure(70);
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,R_all_dB);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[R_all_dB;R_all_dB(1,:)]);
end
hold on
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f70, 'Q7_ACF_all', '-dpng', '-r300')

%% PPI - PSD
f71=figure(71);
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,Pow_all_dB);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[Pow_all_dB;Pow_all_dB(1,:)]);
end
hold on
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f71, 'Q7_PSD_all', '-dpng', '-r300')

%% PPI - Old
f72=figure(72); 
set(gcf,'render','painters');
naz_max = 360;
if length(az_set)<naz_max
	pcolor(x,y,pow_h);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_h;pow_h(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Horizontal Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;
% save
print(f72, 'Q7_power', '-dpng', '-r300')
