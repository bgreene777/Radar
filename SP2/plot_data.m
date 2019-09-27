% 
% Simple program to calculate power and plot for PAR data 
% *** it is assumed that you load the data before running *** 
% 
% Calculate power using the mean of absolute value squared 
% Position indicator (in km)... convert polar (az and range) to rectangular 

%For dual-pol radars, select which polarization to view here
if ~exist('X')
        fprintf('Dual Pol system, changing X_h/X_v to X...\n');
%         X = X_h;
        X = X_v;
end

%The PAR data is in int16 format to save space, this if statement
% will convert the int16 to double format. This is an important
% step as basic matlab commands cannot be performed on int16 data
if ~isfloat(X)
	fprintf('Converting data into floating point numbers ...\n');
	X = double(X);
end
%There has been a little change on the name of the radar OU-PRIME
if strcmp(radar,'OUPRIME'), radar = 'OU-PRIME'; end

el_rad = el/180*pi; 
[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

R0=mean(X(:,:,1:num_pulses).*conj(X(:,:,1:num_pulses)),3);
% R1=mean(X(:,:,2:num_pulses).*conj(X(:,:,1:num_pulses-1)),3);
% R2=mean(X(:,:,3:num_pulses).*conj(X(:,:,1:num_pulses-2)),3);


%%%%%%%%%% Calculate The Power %%%%%%%%%%%%%%%%%%
pow=10*log10(R0);


%%%%%%%%%%%% Plot the Data %%%%%%%%%%%%%%%%%%%%%%
figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
if length(az_set)<naz_max
	pcolor(x,y,pow);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[pow;pow(1,:)]);
end;
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12}',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees']);
colorbar;

%%%%%%%%%%%%  Produce map layover %%%%%%%%%%%%%%%%%%
% bmapover(gca,[],radar,{'OK'});

clear X;


