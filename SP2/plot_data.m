% Brian R. Greene
% Simple program to calculate power and plot for PAR data 
% *** it is assumed that you load the data before running *** 
% 
% Calculate power using the mean of absolute value squared 
% Position indicator (in km)... convert polar (az and range) to rectangular 

close all
clc

%For dual-pol radars, select which polarization to view here
if ~exist('X')
        fprintf('Dual Pol system, changing X_h/X_v to X...\n');
        X = X_h;
%         X = X_v;
end

%The PAR data is in int16 format to save space, this if statement
% will convert the int16 to double format. This is an important
% step as basic matlab commands cannot be performed on int16 data
if ~isfloat(X)
	fprintf('Converting data into floating point numbers ...\n');
	X = double(X);
    X_h = double(X_h);
    X_v = double(X_v);
end
%There has been a little change on the name of the radar OU-PRIME
if strcmp(radar,'OUPRIME'), radar = 'OU-PRIME'; end

el_rad = el/180*pi; 
[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

% R0=mean(X(:,:,1:num_pulses).*conj(X(:,:,1:num_pulses)),3);
R_h=mean(X_h(:,:,1:num_pulses).*conj(X_h(:,:,1:num_pulses)),3);
R_v=mean(X_v(:,:,1:num_pulses).*conj(X_v(:,:,1:num_pulses)),3);


%%%%%%%%%% Calculate The Power %%%%%%%%%%%%%%%%%%
% pow=10*log10(R0);
pow_h = 10*log10(R_h);
pow_v = 10*log10(R_v);

% remove low power returs
i_low_h = find(pow_h < 20);
i_low_v = find(pow_v < 20);

pow_h(i_low_h) = nan;
pow_v(i_low_v) = nan;

% diff
pow_diff = pow_h - pow_v;

%%%%%%%%%%%% Plot the Data %%%%%%%%%%%%%%%%%%%%%%
% Horizontal polarization
figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
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

% Vertical polarization
figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
if length(az_set)<naz_max
	pcolor(x,y,pow_v);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_v;pow_v(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Vertical Power ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;

% Difference
figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
if length(az_set)<naz_max
	pcolor(x,y,pow_diff);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow_diff;pow_diff(1,:)]);
end
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12} Difference between h and v ',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;

%%%%%%%%%%%%  Produce map layover %%%%%%%%%%%%%%%%%%
% bmapover(gca,[],radar,{'OK'});

clear X;


