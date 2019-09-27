% BMAPOVER   Boon Leng's Map overlay
%    BMAPOVER draws a map overlay on top of the current axis, assuming the
%    the KOUN radar (97.463161W,35.236208N) as the origin (0,0).
%    The axis limit is in units of km.
%
%    More usage options:
%    BMAPOVER(AXES_H) draws map overlay on AXES_H
%
%    BMAPOVER(AXES_H,FLAGS,ORIGIN,ONESTATE) draws the map on AXES_H (can 
%    be a vector) and uses FLAGS to determine which lines to draw.
%    FLAGS is a logical array recognized as:
%    FLAGS = [State, County, Interstate, Local HWY, Capital, County Seat,
%             Big Town, Small Town, Others, OK Mesonet, Weather Station, GUI]
%    SUBSTATE = 0 or 1 to plot only the viewing state only.
%    SUBSTATE = 'OK' to plot only OK state
%    SUBSTATE = {'OK','KS'} to plot OK and KS states
%
%    Usage: bmapover(gca) draws using default options; or
%           bmapover(gca,[1 0 1]) draws state and interstate only
%           bmapover(gca,[1 3]) is the same as above
%           bmapover(gca,[],'KCYR') draws with origin as 'KCYR'
%           bmapover(gca,[],'KDDC',{'NE','KS'}) draws with just NE and KS states
%           bmapover(gca,[],{-104,36,'Radar','OK'}) plots around (lon,lat)
%
%    VIS = BMAPOVER returns a structure containing the mapping details
%
%
%    Created on 6/5/2005
%    Last updated on 
%    Boon Leng Cheong
%
%    Version 1.00 - Updated on 6/17/2010
%                 - Added full WSR-88D network.
%    Version 0.95 - Updated on 5/18/2008
%                 - Added KFDR.
%    Version 0.95 - Updated on 5/15/2006
%                 - UI is independent from desktop variable VIS.
%                 - New map data from http://www.nationalatlas.gov
%                 - Support 48 continental states for state and county
%                   borders, state capital and interstate highway
%                 - New map data includes Mesonet stations
%                 - Airport labels are removed
%                 - Labels are not drawn first so less objects to handle
%    Version 0.94 - Took out the order option, not used anyway.
%                 - Added flags option, don't plot objectss that are not
%                   requested.  Sometimes this speed things up by 4x !
%    Version 0.93 - Improve compatibility with Matlab v6
%    Version 0.92 - No points for city labels.
%    Version 0.91 - Use rotations, rotate (-longitude) along z-axis, 
%                   then (-latitude) along y-axis, look along x-axis
%                   toward the origin. z = north, y towards greenwich 
%    Version 0.9  - More flexible, uses own axes
%                 - Supports multiple axes, AXES_H can be a vector
%    Version 0.1  - Buggy, use at your own risk
%
function [VIS] = bmapover(axes_handler_set,flags,origin,substate)

if ~exist('axes_handler_set','var') || isempty(axes_handler_set)
    if isempty(get(gcf,'CurrentAxes')), 
        axes_handler_set = axes('Units','Normalized','Position',[0.08 0.06 0.9 0.9]);
        axis(0.5*[-2350 2450 -1500 2100]);
    else
        axes_handler_set = gca;
    end
end
if exist('flags','var') && ~isempty(flags) && any(flags>1),
    tmp = logical([zeros(1,11) 1]); tmp(flags) = 1; flags = logical(tmp);
end
if ~exist('flags','var') || isempty(flags), flags = logical([ones(1,10) 1 1]); end
if (length(flags)~=12), flags(length(flags)+1:11) = 0; flags(12) = 1; end
if ~exist('origin','var') || isempty(origin), origin = 'OU-PRIME'; end
if ~exist('substate','var') || isempty(flags), substate = 0; end
if ischar(substate) || iscell(substate)
    tmp = {'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DE', 'FL', 'GA', 'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA',...
           'MA', 'MD', 'ME', 'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY',...
           'OH', 'OK', 'OR', 'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY'};
    if any(~ismember(substate,tmp))
        fprintf('One of the specified states is not available')
        return
    end
    ref_st = substate;
    substate = true;
end

% load('bmapover_map_full.mat');
load('bmapover_map.mat');

r_earth = 6378.1;  % The Earth's Radius = 6378.1 km
% Reminder: When adding a station, lab, lon, lat and st must all be filled
% 'CHILL' -104.63708  40.44625   CO
%{'Blue Canyon', -98.54752, 34.85820, 'OK'};

%     origin = {-97.930556,34.81278,'KRSP','OK'};
%     origin = {-97.956111,35.03139,'KSAO','OK'};

% MAP(11).lab  = {'KOUN/PAR','KDDC','KPUX','KTLX','KFTG','XPOL','KCYR','KRSP','KSAO',...
% 	'KLWE','KFDR','KDYX','KTFX','OU-PRIME','KCRP','KBUF','KTYX'};
% MAP(11).lon = [-97.463161 -99.968800 -104.181638 -97.277831 -104.545791 -104.63824 -98.25212 -97.93056 -97.95611 ...
% 	-98.27201 -98.97611 -99.25417 -111.38444 -97.433701 -97.51083 -78.73694 -75.68000];
% MAP(11).lat = [ 35.236208  37.760804   38.459453  35.333411   39.786593   40.44616  34.87398  34.81278  35.03139 ...
% 	 34.62381  34.36222  32.53833   47.45972  35.180099  27.78389  42.94861  43.75583];
% MAP(11).st = {'OK','KS','CO','OK','CO','CO','OK','OK','OK',...
% 	'OK','OK','TX','MT','OK','TX','NH','NY'};

MAP(11).lab = {'OU-PRIME','NWRT','KCYR','KRSP','KSAO', 'PX-1000', 'KLWE',...
	'KABR'     ,'KENX'     ,'KABX'     ,'KFDR'     ,'KAMA'     ,'PAHG'     ,'PGUA'     ,'KFFC'     ,'KEWX'     ,'KBBX'     ,'PABC'     ,'KBLX'     ,'KBGM'     ,'KBMX'     ,'KBIS'     ,'KCBX'     ,'KBOX'     ,'KBRO'     ,'KBUF'     ,'KCXX'     ,'RKSG'     ,'KFDX'     ,'KICX'     ,'KCLX'     ,'KRLX'     ,'KCYS'     ,'KLOT'     ,'KILN'     ,'KCLE'     ,'KCAE'     ,'KGWX'     ,'KCRP'     ,'KFWS'     ,'KDVN'     ,'KFTG'     ,'KDMX'     ,'KDTX'     ,'KDDC'     ,'KDOX'     ,'KDLH'     ,'KDYX'     ,'KEYX'     ,'KEVX'     ,'KEPZ'     ,'KLRX'     ,'KBHX'     ,'PAPD'     ,'KFSX'     ,'KHPX'     ,'KGRK'     ,'KPOE'     ,'KEOX'     ,'KSRX'     ,'KIWX'     ,'KAPX'     ,'KGGW'     ,'KGLD'     ,'KMVX'     ,'KGJX'     ,'KGRR'     ,'KTFX'     ,'KGRB'     ,'KGSP'     ,'KRMX'     ,'KUEX'     ,'KHDX'     ,'KCBW'     ,'KHGX'     ,'KHTX'     ,'KIND'     ,'KJKL'     ,'KJAN'     ,'KJAX'     ,'RODN'     ,'PHKN'     ,'KEAX'     ,'KBYX'     ,'PAKC'     ,'KMRX'     ,'RKJK'     ,'KARX'     ,'LPLA'     ,'KLCH'     ,'KESX'     ,'KDFX'     ,'KILX'     ,'KLZK'     ,'KVTX'     ,'KLVX'     ,'KLBB'     ,'KMQT'     ,'KMXX'     ,'KMAX'     ,'KMLB'     ,'KNQA'     ,'KAMX'     ,'PAIH'     ,'KMAF'     ,'KMKX'     ,'KMPX'     ,'KMBX'     ,'KMSX'     ,'KMOB'     ,'PHMO'     ,'KVAX'     ,'KMHX'     ,'KOHX'     ,'KLIX'     ,'KOKX'     ,'PAEC'     ,'KAKQ'     ,'KLNX'     ,'KTLX'     ,'KOAX'     ,'KPAH'     ,'KPDT'     ,'KDIX'     ,'KIWA'     ,'KPBZ'     ,'KSFX'     ,'KGYX'     ,'KRTX'     ,'KPUX'     ,'KRAX'     ,'KUDX'     ,'KRGX'     ,'KRIW'     ,'KFCX'     ,'KJGX'     ,'KDAX'     ,'KLSX'     ,'KMTX'     ,'KSJT'     ,'KNKX'     ,'KMUX'     ,'KHNX'     ,'TJUA'     ,'KSOX'     ,'KATX'     ,'KSHV'     ,'KFSD'     ,'PACG'     ,'PHKI'     ,'PHWA'     ,'KOTX'     ,'KSGF'     ,'KCCX'     ,'KLWX'     ,'KTLH'     ,'KTBW'     ,'KTWX'     ,'KEMX'     ,'KINX'     ,'KVNX'     ,'KVBX'     ,'KICT'     ,'KLTX'     ,'KYUX'     ,...
    'KTYX','RaXpol_hurricane'};
MAP(11).lon = [-97.43550, -97.463161, -98.25212 -97.93056 -97.95611, -98.27201,...
	-98.413056, -74.063889,-106.823889, -98.976389,-101.709167,-151.351389,-144.811389, -84.565833, -98.028333,-121.631667,-161.876389,-108.606667, -75.984722, -86.769722,-100.760556,-116.235556, -71.136944, -97.418889, -78.736667, -73.166944,-127.021111,-103.630000,-112.862222, -81.042222, -81.723056,-104.806111, -88.084722, -83.821667, -81.859722, -81.118333, -88.328889, -97.511111, -97.303056, -90.580833,-104.545833, -93.722778, -83.471667, -99.968889, -75.439722, -92.209722, -99.254444,-117.560833, -85.921389,-106.698056,-116.802778,-124.291944,-147.501667,-111.197778, -87.285000, -97.383056, -92.975833, -85.459444, -94.361667, -85.700000, -84.719722,-106.625000,-101.700278, -97.325556,-108.213889, -85.544722,-111.385278, -88.111389, -82.220000, -75.457778, -98.441944,-106.122778, -67.806389, -95.079167, -86.083333, -86.280278, -83.313056, -90.080000, -81.701944,-127.909722,-155.777778, -94.264444, -81.703056,-156.629444, -83.401667,-126.622222, -91.191111, -27.321667, -93.215833,-114.891389,-100.280556, -89.336944, -92.262222,-119.179444, -85.943889,-101.814167, -87.548333, -85.789722,-122.717222, -80.654167, -89.873333, -80.412778,-146.303056,-102.189167, -88.550556, -93.565556,-100.865000,-113.986111, -88.239722,-157.180000, -83.001667, -76.876111, -86.562500, -89.825556, -72.863889,-165.295000, -77.007222,-100.576389, -97.277778, -96.366667, -88.771944,-118.852778, -74.410833,-111.670000, -80.218056,-112.686111, -70.256389,-122.965278,-104.181389, -78.489722,-102.829722,-119.462222,-108.477222, -80.273889, -83.351111,-121.677778, -90.682778,-112.447778,-100.492500,-117.041944,-121.898333,-119.632222, -66.078056,-117.635833,-122.495833, -93.841389, -96.729444,-135.529167,-159.552222,-155.568889,-117.626667, -93.400556, -78.003611, -77.477778, -84.328889, -82.401667, -96.232500,-110.630278, -95.564722, -98.127778,-120.396944, -97.442778, -78.428889,-114.656667,...
	-75.680000,-76.6610];
MAP(11).lat = [35.180299, 35.236208, 34.87398  34.81278  35.03139,  34.62381,...
     45.455833,  42.586389,  35.149722,  34.362222,  35.233333,  60.725833,  13.452500,  33.363611,  29.703889,  39.496111,  60.791944,  45.853889,  42.199722,  33.172222,  46.770833,  43.490556,  41.955833,  25.916111,  42.948889,  44.511111,  36.955833,  34.635278,  37.590833,  32.655556,  38.311111,  41.151944,  41.604722,  39.420278,  41.413056,  33.948611,  33.896667,  27.784167,  32.573056,  41.611667,  39.786667,  41.731111,  42.699722,  37.760833,  38.825556,  46.836944,  32.538333,  35.097778,  30.564444,  31.873056,  40.739722,  40.498333,  65.035000,  34.574444,  36.736667,  30.721944,  31.155556,  31.460556,  35.290556,  41.358889,  44.907222,  48.206389,  39.366944,  47.527778,  39.062222,  42.893889,  47.459722,  44.498333,  34.883333,  43.467778,  40.320833,  33.076389,  46.039167,  29.471944,  34.930556,  39.707500,  37.590833,  32.317778,  30.484722,  26.301944,  20.125556,  38.810278,  24.597500,  58.679444,  36.168611,  35.924167,  43.822778,  38.730278,  30.125278,  35.701111,  29.272778,  40.150556,  34.836389,  34.411667,  37.975278,  33.654167,  46.531111,  32.536667,  42.081111,  28.113333,  35.344722,  25.611111,  59.461389,  31.943333,  42.967778,  44.848889,  48.392500,  47.041111,  30.679444,  21.132778,  30.390278,  34.776111,  36.247222,  30.336667,  40.865556,  64.511389,  36.983889,  41.957778,  35.333056,  41.320278,  37.068333,  45.690556,  39.946944,  33.289167,  40.531667,  43.105833,  43.891389,  45.714722,  38.459444,  35.665556,  44.125000,  39.754167,  43.066111,  37.024444,  32.675278,  38.501111,  38.698889,  41.262778,  31.371389,  32.918889,  37.155278,  36.314167,  18.115556,  33.817778,  48.194444,  32.450833,  43.587778,  56.852778,  21.894167,  19.095000,  47.680278,  37.235278,  40.923056,  38.975278,  30.397500,  27.705556,  38.996944,  31.893611,  36.175000,  36.740833,  34.838056,  37.654722,  33.989444,  32.495278,...
	 43.75583,34.76330];
MAP(11).st = {'OK','OK','OK','OK','OK','OK',...
	'SD'       ,'NY'       ,'NM'       ,'OK'       ,'TX'       ,'AK'       ,'GU'       ,'GA'       ,'TX'       ,'CA'       ,'AK'       ,'MT'       ,'NY'       ,'AL'       ,'ND'       ,'ID'       ,'MA'       ,'TX'       ,'NY'       ,'VT'       ,'ea'       ,'NM'       ,'UT'       ,'SC'       ,'WV'       ,'WY'       ,'IL'       ,'OH'       ,'OH'       ,'SC'       ,'MS'       ,'TX'       ,'TX'       ,'IA'       ,'CO'       ,'IA'       ,'MI'       ,'KS'       ,'DE'       ,'MN'       ,'TX'       ,'CA'       ,'FL'       ,'TX'       ,'NV'       ,'CA'       ,'AK'       ,'AZ'       ,'KY'       ,'TX'       ,'LA'       ,'AL'       ,'AR'       ,'IN'       ,'MI'       ,'MT'       ,'KS'       ,'ND'       ,'Co'       ,'MI'       ,'MT'       ,'WI'       ,'SC'       ,'NY'       ,'NE'       ,'NM'       ,'ME'       ,'TX'       ,'AL'       ,'IN'       ,'KY'       ,'MS'       ,'FL'       ,'wa'       ,'HI'       ,'MO'       ,'FL'       ,'AK'       ,'TN'       ,'ea'       ,'WI'       ,'es'       ,'LA'       ,'NV'       ,'TX'       ,'IL'       ,'AR'       ,'CA'       ,'KY'       ,'TX'       ,'MI'       ,'AL'       ,'OR'       ,'FL'       ,'TN'       ,'FL'       ,'AK'       ,'TX'       ,'WI'       ,'MN'       ,'ND'       ,'MT'       ,'AL'       ,'HI'       ,'GA'       ,'NC'       ,'TN'       ,'LA'       ,'NY'       ,'AK'       ,'VA'       ,'NE'       ,'OK'       ,'NE'       ,'KY'       ,'OR'       ,'PA'       ,'AZ'       ,'PA'       ,'ID'       ,'ME'       ,'OR'       ,'CO'       ,'NC'       ,'SD'       ,'NV'       ,'WY'       ,'VA'       ,'GA'       ,'CA'       ,'MO'       ,'UT'       ,'TX'       ,'CA'       ,'CA'       ,'CA'       ,'PR'       ,'CA'       ,'WA'       ,'LA'       ,'SD'       ,'AK'       ,'HI'       ,'HI'       ,'WA'       ,'MO'       ,'PA'       ,'VA'       ,'FL'       ,'FL'       ,'KS'       ,'AZ'       ,'OK'       ,'OK'       ,'CA'       ,'KS'       ,'NC'       ,'AZ'       ,...
	'NY','NC'};

if iscell(origin)
    if numel(origin)==4
        MAP(11).lab{end+1} = origin{3};
        MAP(11).lon(end+1) = origin{1};
        MAP(11).lat(end+1) = origin{2};
        MAP(11).st{end+1} = origin{4};
        origin = origin{3};
    else
        fprintf('Invalid input structure, should be {lat,lon,''Label'',''ST''}.\n');
        return;
    end
end
if ismember(origin,{'KCRI','PAR','KOUN'}), origin = 'NWRT'; end
if ~ischar(origin)&&~isstruct(origin), error('Origin must be either a string or a structure'); end
if ~ismember(origin,MAP(11).lab), error(['No origin: ',origin]); end

[yesno,pos] = ismember(origin,MAP(11).lab); %#ok<ASGLU>
lon_o = MAP(11).lon(pos);
lat_o = MAP(11).lat(pos);

% Get the properties of the plotting axes
axes_handler = axes_handler_set(1);
axes(axes_handler) %#ok<*MAXES>
alim = [get(axes_handler,'XLim'),get(axes_handler,'YLim')];  % Axis limit in m
% Initialize some handles before plotting the map
num_ax = length(axes_handler_set);
tmp = -1*ones(num_ax,6);
VIS = struct('fh',[],'peer',[],'ax',axes_handler_set,'h',tmp,'ht1',tmp(:,1),'ht2',tmp(:,1),'ht3',tmp(:,1),...
             'ht4',tmp(:,1),'ht5',tmp(:,1),'ht6',tmp(:,1),'ht7',tmp(:,1));
% Coverage in km to plot, the buffer is a bit tricky, some segments are 50km long!
if max(alim(2)-alim(1),alim(4)-alim(3))>100
    tmp = alim+0.3*(alim(2)-alim(1))*[-1 1 -1 1];
else
    tmp = alim+40*[-1 1 -1 1];
end
% Another extra a little bit to work with
LL = [tmp(1:2)/cos(lat_o/180*pi) tmp(3:4)]./r_earth*180/pi+[lon_o lon_o lat_o lat_o]+...
     [max(0.1*(tmp(2)-tmp(1))/cos(lat_o/180*pi),10)*[-1 1] max(0.1*(tmp(4)-tmp(3)),10)*[-1 1]]./r_earth*180/pi;
% Pre-processing to focus on needed subsets
if ((LL(1)>-118)||(LL(2)<-74)||(LL(3)>30)||(LL(4)<45))
    lon_v = 0.5*(LL(1)+LL(2)); lat_v = 0.5*(LL(3)+LL(4));
    if ~exist('ref_st','var')
        % Find a state bounding box that encloses the viewing center and use it as ref_st
        yn = (lon_v>MAP(2).box(1,:))&(lon_v<MAP(2).box(2,:))&(lat_v>MAP(2).box(3,:))&(lat_v<MAP(2).box(4,:));
        % idx = find(yn); ref_st = unique(MAP(2).st(idx));
		ref_st = unique(MAP(2).st(yn));
        %if (length(ref_st)>1), fprintf('I''m guessing the reference is %s\n',ref_st{1}); end
        if isempty(ref_st), fprintf('Sorry out of data.\n'); return; end
    end
    if (~any(ismember(ref_st,{'OK'}))&&substate), flags(10) = 0; end
    
    for idx = find(flags(1:4))
        if substate
            yn = ismember(MAP(idx).st,ref_st);
        else
            yn = isoverlap(LL,MAP(idx).box);
        end
        MAP(idx).lon = [MAP(idx).lon{yn}];
        MAP(idx).lat = [MAP(idx).lat{yn}];
        MAP(idx).st = MAP(idx).st(yn);
    end
    if substate
        for idx = 4+find(flags(5:11))
            yn = ismember(MAP(idx).st,ref_st);
            MAP(idx).lon = MAP(idx).lon(yn);
            MAP(idx).lat = MAP(idx).lat(yn);
            MAP(idx).lab = MAP(idx).lab(yn);
            MAP(idx).st = MAP(idx).st(yn);
        end
    end
else
    % Convert cell to double array
    for idx = 1:4
        MAP(idx).lon = [MAP(idx).lon{:}]; MAP(idx).lat = [MAP(idx).lat{:}];
    end
end

% Convert longitude/latitude to 3-D on the earth surface as a sphere, then
% rotate along z-axis, then along y-axis.  View the map from x-axis, which 
% is viewing Y vs. Z as East vs. North.
%thx = 0/180*pi;
thy = lat_o/180*pi;
thz = -lon_o/180*pi;
%Rx = [1 0 0; 0 cos(thx) -sin(thx); 0 sin(thx) cos(thx)];
Ry = [cos(thy) 0 sin(thy); 0 1 0; -sin(thy) 0 cos(thy)];
Rz = [cos(thz) -sin(thz) 0; sin(thz) cos(thz) 0; 0 0 1];
R = Ry*Rz;

% Earth's equitorial radius a = 6,378.137 km
% a = 6378.137;
% Earth's polar radius b = 6,356.7523 km
% b = 6356.7523;

% r_earth = sqrt(((a^2+cos(phi)).^2+(b^2*sin(phi)).^2)./((a*cos(phi)).^2+(b*sin(phi)).^2));
for idx = find(flags(1:4))
	
    % Absolute position for the points/lines
    xyz = r_earth*[cos(MAP(idx).lat*pi/180).*cos(MAP(idx).lon*pi/180); ...
                   cos(MAP(idx).lat*pi/180).*sin(MAP(idx).lon*pi/180); ...
                   sin(MAP(idx).lat*pi/180)];
% 	phi = MAP(idx).lat*pi/180;
% 	theta = MAP(idx).lon*pi/180;
% 	r_earth = sqrt(((a^2+cos(phi)).^2+(b^2*sin(phi)).^2)./((a*cos(phi)).^2+(b*sin(phi)).^2));
% 	xyz = ([1 1 1]'*r_earth).*[cos(phi).*cos(theta); cos(phi).*sin(theta); sin(phi)];

    if size(xyz,1)~=3
        MAP(idx).x = [];
        MAP(idx).y = [];
        %MAP(idx).z = [];
    else
        MAP(idx).x = R(2,:)*xyz;
        MAP(idx).y = R(3,:)*xyz;
        %MAP(idx).z = R(1,:)*xyz-r_earth;
    end
    % Within plotting domain or NaN (need them for different lines)
    mask = (MAP(idx).x>=tmp(1)&MAP(idx).x<=tmp(2)&...
            MAP(idx).y>=tmp(3)&MAP(idx).y<=tmp(4))|...
           ~isfinite(MAP(idx).x);
    mask2 = ~isfinite(MAP(idx).x);
    stidx = [1 find(mask2)];
    for jdx = 1:length(stidx)-1
        % If only some segments of the line are inside
        if any(mask(stidx(jdx)+1:stidx(jdx+1)-1)) && ~all(mask(stidx(jdx)+1:stidx(jdx+1)-1)),
            % Points between NaN that are outside plotting domain
            eidx = find(~mask(stidx(jdx)+1:stidx(jdx+1)-1));
            mask(stidx(jdx)+eidx) = 1;
            MAP(idx).x(stidx(jdx)+eidx) = nan;
        end
    end
    if all(~mask)
        % Put a dummy point if there is nothing to plot
        MAP(idx).x = nan;
        MAP(idx).y = nan;
        %MAP(idx).z = nan;
        flags(idx) = 0;
    else
        % Extract the okay points
        MAP(idx).x = MAP(idx).x(mask);
        MAP(idx).y = MAP(idx).y(mask);
        %MAP(idx).z = MAP(idx).z(mask);
    end
    % Take out consecutive NaN
    mask = [true,isfinite(MAP(idx).x(1:end-1))|isfinite(MAP(idx).x(2:end))];
    MAP(idx).x = MAP(idx).x(mask);
    MAP(idx).y = MAP(idx).y(mask);
    %MAP(idx).z = MAP(idx).z(mask);
end

% For labels: Extract the labels that might be used later
for idx = 4+find(flags(5:11))

    xyz = r_earth*[cos(MAP(idx).lat*pi/180).*cos(MAP(idx).lon*pi/180); ...
                   cos(MAP(idx).lat*pi/180).*sin(MAP(idx).lon*pi/180); ...
                   sin(MAP(idx).lat*pi/180)];
% 	phi = MAP(idx).lat*pi/180;
% 	theta = MAP(idx).lon*pi/180;
% 	r_earth = sqrt(((a^2+cos(phi)).^2+(b^2*sin(phi)).^2)./((a*cos(phi)).^2+(b*sin(phi)).^2));
% 	xyz = ([1 1 1]'*r_earth).*[cos(phi).*cos(theta); cos(phi).*sin(theta); sin(phi)];

    x = R(2,:)*xyz;
    y = R(3,:)*xyz;
    %z = R(1,:)*xyz-r_earth;
    mask = (x>=tmp(1))&(x<=tmp(2))&(y>=tmp(3))&(y<=tmp(4));
    % Special case: copy all the mesonet stations out
    if (idx==10), VIS.mnet_x = x; VIS.mnet_y = y; VIS.mnet_stid = MAP(idx).lab; end
    if isempty(mask) || all(~mask)
        % Put a dummy point if there is nothing to plot
        MAP(idx).x = nan;
        MAP(idx).y = nan;
        %MAP(idx).z = nan;
        MAP(idx).lab = {''};
        MAP(idx).st = {''};
        flags(idx) = 0;
    else
        MAP(idx).x = x(mask);
        MAP(idx).y = y(mask);
        %MAP(idx).z = z(mask)-r_earth;
        MAP(idx).lab = MAP(idx).lab(mask);
        MAP(idx).st = MAP(idx).st(mask);
    end
end

br = get(gca,'Color');
br = [0.3 0.59 0.11]*br(:);      % Brightness of gca's color

if br>0.5
    clr = [0.7843  0.7451  0.4706; ... % State
           0.4000  0.4000  0.4000; ... % County
           0.6275  0.5490  0.4706; ... % Interstate
           0.9804  0.7059  0.2118; ... % Local HWY
                0       0       0; ... % State Capital
           0.6000  0.3000  1.0000; ... % KOUN/PAR, etc
           0.3500  0.2800  0.2000; ... % County Seat
           0.3000  0.2000  0.1000; ... % Pop>10k
           0.3000  0.3000  0.3000; ... % Pop>1000
           0.3000  0.3000  0.3000; ... % Others
           0.3500  0.2000       0];    % Mesonet
else
    clr = [0.9216  0.6549       0; ... % State
           0.6627  0.3382  0.1922; ... % County
           1.0000       0       0; ... % Interstate
           0.5000  0.1000  0.1000; ... % Local HWY
           1.0000  1.0000  1.0000; ... % State Capital
           0.7300  0.5250  1.0000; ... % KOUN/PAR, etc
           0.7843  0.7843  0.7843; ... % County Seat
           0.6275  0.6275  0.6275; ... % Pop>10k
           0.4706  0.4706  0.4706; ... % Pop>1000
           0.3137  0.3137  0.3137; ... % Others
           1.0000  0.9075  0.8151];    % Mesonet
end

% Drawing
for iax=1:num_ax
    axes(axes_handler_set(iax)) %#ok<*LAXES>
    hold on
    if flags(2), VIS.h(iax,2) = plot(MAP(2).x,MAP(2).y,'Color',clr(2,:)); end                 % County
    if flags(4), VIS.h(iax,4) = plot(MAP(4).x,MAP(4).y,'Color',clr(4,:),'LineWidth',1); end   % Local HWY
    if flags(3), VIS.h(iax,3) = plot(MAP(3).x,MAP(3).y,'Color',clr(3,:),'LineWidth',2); end   % Interstate
    if flags(1), VIS.h(iax,1) = plot(MAP(1).x,MAP(1).y,'Color',clr(1,:),'LineWidth',2); end   % State
    if flags(10), VIS.h(iax,5:4+double(~isempty(MAP(10).x))) = plot(MAP(10).x,MAP(10).y,'^',...
        'Clipping','On','Color',clr(11,:),'MarkerFaceColor',clr(11,:),'MarkerSize',4); end    % OK Mesonet
    if flags(5), VIS.ht2(iax,1:length(MAP(5).x)) = text(MAP(5).x,MAP(5).y,MAP(5).lab,...
        'Color',clr(5,:),'FontWeight','Bold','FontSize',12,'Clipping','On',...
        'HorizontalAlignment','Center','VerticalAlignment','Middle');
    end                     % State Capital
    if flags(11), VIS.ht1(iax,1:length(MAP(11).lab)) = text(MAP(11).x,MAP(11).y,MAP(11).lab,...
        'Color',clr(6,:),'FontWeight','Bold','FontSize',10,'Clipping','On',...
        'HorizontalAlignment','Left','VerticalAlignment','Bottom');                           % KOUN/PAR, etc
    end
    if flags(11), VIS.h(iax,6) = plot(MAP(11).x,MAP(11).y,'o',...
        'Color',clr(6,:),'Markersize',3,'MarkerFaceColor',clr(6,:),'Clipping','On');          % KOUN/PAR, etc
    end
    VIS.ht3(iax,1) = text(0,0,'','Color',clr(7,:));
    VIS.ht4(iax,1) = text(0,0,'','Color',clr(8,:));
    VIS.ht5(iax,1) = text(0,0,'','Color',clr(9,:));
    VIS.ht6(iax,1) = text(0,0,'','Color',clr(10,:));
    VIS.ht7(iax,1) = text(0,0,'','Color',clr(11,:));
    hold off
    set(gca,'XLim',alim(1:2),'YLim',alim(3:4),'DataAspect',[1 1 1],'Layer','Top')
end
% Decide zoom level if no flags were supplied, hide some lines. Labels are a bit tricky at this point
%if (nargin<2)
    scl = [250 100];
    cvg = 0.5*(alim(2)-alim(1))+0.5*(alim(4)-alim(3));
    if (cvg>scl(1))  % Almost entire state
        if flags(4), set(VIS.h(:,4),'Visible','Off'); end
        if flags(10), set(VIS.h(:,5),'Visible','Off'); end
    end
%end

% A GUI for showing/hiding plotted objects
if flags(12)
    if br>0.5, bgcolor = '[1 1 1]'; else bgcolor = ['[ ',num2str(get(gca,'color'),'%.3f '),']']; end
    VIS.peer = get(axes_handler_set(1),'Parent');
    VIS.fh = VIS.peer.Number+99;
    tmp = get(VIS.peer,'DeleteFcn');
    if isempty(tmp)
        set(VIS.peer,'DeleteFcn',['if ishandle(',num2str(VIS.fh),'); delete(',num2str(VIS.fh),'); end']);
    else
        tmp = [tmp '; if ishandle(',num2str(VIS.fh),'); delete(',num2str(VIS.fh),'); end'];
        set(VIS.peer,'DeleteFcn',tmp);
    end
    figure(VIS.fh)
    clf
    ibn = 1; uiysize = sum(flags([1:10 10 11]))*25+15;
    set(gcf,'Menubar','None','NumberTitle','Off','Color',[0.842 0.794 0.721]) % [0.65 0.55 0.45]
    tmp = get(gcf,'Position'); set(gcf,'Position',[tmp(1:2) 150 uiysize+5]);
    figui = get(VIS.fh,'Position'); figui = figui([1 3 2 4])+[0 figui(1) 0 figui(2)];
    figax = get(VIS.peer,'Position'); figax = figax([1 3 2 4])+[0 figax(1) 0 figax(2)];
    if isoverlap(figax,figui), set(gcf,'Position',[figax(1)-155 figui(3) 150 uiysize+5]); end
    for idx = find(flags([1:10 10 11]))
        switch idx
            case 1, lab = 'State Borders'; hdl = 'VIS.h(:,1)'; fs = '11'; w = 'Bold';
            case 2, lab = 'County Borders'; hdl = 'VIS.h(:,2)'; fs = '11'; w = 'Bold';
            case 3, lab = 'Interstate'; hdl = 'VIS.h(:,3)'; fs = '11'; w = 'Bold';
            case 4, lab = 'Local HWY'; hdl = 'VIS.h(:,4)'; fs = '11'; w = 'Bold';
            case 5, lab = 'State Capital'; hdl = 'VIS.ht2'; fs = '11'; w = 'Bold';
            case 6, lab = 'County Seat'; hdl = 'VIS.ht3'; fs = '10'; w = 'Normal';
            case 7, lab = 'Pop>10k'; hdl = 'VIS.ht4'; fs = '10'; w = 'Normal';
            case 8, lab = 'Pop>1000'; hdl = 'VIS.ht5'; fs = '9'; w = 'Normal';
            case 9, lab = 'Others'; hdl = 'VIS.ht6'; fs = '9'; w = 'Normal';
            case 10, lab = 'OK MESONET'; hdl = 'VIS.h(:,5)'; fs = '7'; w = 'Normal';
            case 11, lab = 'OK MESONET ID'; hdl = 'VIS.ht7'; fs = '7'; w = 'Normal';
            case 12, lab = 'Radar'; hdl = '[VIS.ht1(:); VIS.h(:,6)]'; fs = '10'; w = 'Bold';
        end
        if (idx<=10), lset = num2str(idx); else lset = num2str(idx-1); end
        eval(['c = ',hdl,'; if ishandle(c(1)), c = sprintf(''[%.3f %.3f %.3f]'',0.9*get(c(1),''Color'')); ',...
              'else, c = ''[0 0 0]'', end']);
        comm = ['VIS.hb(',num2str(ibn),') = uicontrol(''Style'',''PushButton'',''Unit'',''Pixel'',',...
                '''Position'',[15 ',num2str(uiysize-25*ibn),' 120 20],''String'',''',lab,''',',...
                '''FontSize'',',fs,',''ForegroundColor'',',c,',''FontWeight'',''',w,''',',...
                '''BackgroundColor'',',bgcolor,',''Visible'',''On'');'];
        eval(comm);
        if ((idx>=5)&&(idx<=9)) || (idx==11),
            % For text objects, only draw them when it's needed because they take one handle per string!
            if (idx==11), pad1 = 'smll = 0.015*(alims(4)-alims(3)); '; pad2 = '-smll'; else pad1 = ''; pad2 = ''; end
            comm = ['set(VIS.hb(',num2str(ibn),'),'...
            '''Callback'',', ...
            '''a = exist(''''VIS'''',''''var''''); ',...
            'VIS = get(gcf,''''UserData''''); ',...
            'if any(ishandle(',hdl,')), ',...
                'delete(',hdl,'(ishandle(',hdl,'))), ',...
                hdl,' = repmat(1.001,[length(VIS.ax) 1]); ',...
            'else, ',...
              'alims = [get(VIS.ax(1),''''XLim'''') get(VIS.ax(1),''''YLim'''')]; ',...
              'alims = alims+0.05*(alims(2)-alims(1))*[1 -1 1 -1]; ',pad1,...
              'loc = (VIS.MAP(',lset,').x>alims(1))&(VIS.MAP(',lset,').x<alims(2))&',...
                    '(VIS.MAP(',lset,').y>alims(3))&(VIS.MAP(',lset,').y<alims(4)); ',...
              'for iax = 1:length(VIS.ax), ',...
                'axes(VIS.ax(iax)), ',...
                 hdl,'(iax,1:sum(loc)) = ',...
                    'text(VIS.MAP(',lset,').x(loc),VIS.MAP(',lset,').y(loc)',pad2,',VIS.MAP(',lset,').lab(loc),',...
                          '''''Color'''',get(VIS.hb(',num2str(ibn),'),''''ForegroundColor''''),'...
                          '''''FontSize'''',get(VIS.hb(',num2str(ibn),'),''''FontSize''''),',...
                          '''''FontWeight'''',get(VIS.hb(',num2str(ibn),'),''''FontWeight''''),'...
                          '''''Clipping'''',''''On'''', ', ...
                          '''''HorizontalAlignment'''',''''Center'''', ', ...
                          '''''VerticalAlignment'''',''''Middle''''); ', ...
              'end; ',...
            'end; ',...
            'set(VIS.fh,''''UserData'''',VIS); if ~a, clear VIS; end; clear a alims smll iax loc;'');'];
        else
            % For lines, just hide/show them using a toggle switch.
            comm = ['set(VIS.hb(',num2str(ibn),'),'...
                '''Callback'',', ...
                '''VIS = get(gcf,''''UserData''''); ',...
                'if any(strcmp(lower(get(',hdl,',''''Visible'''')),''''on'''')),',...
                'set(',hdl,',''''Visible'''',''''Off''''); ', ...
                'else, set(',hdl,',''''Visible'''',''''On''''); end; clear VIS;'');'];
        end
        eval(comm);
        ibn = ibn+1;
    end
    % Final touch, avoid returning handle = 0 (root handle)
    VIS.hb = VIS.hb(VIS.hb~=0);
    tmp = {'lon','lat','box'};
    VIS.MAP = rmfield(MAP,intersect(tmp,fieldnames(MAP)));
    set(gcf,'UserData',VIS,'Tag','bmapover');
    figure(VIS.peer)
end

% Delete empty text labels, they were used to initialize a process only
delete([VIS.ht3 VIS.ht4 VIS.ht5 VIS.ht6 VIS.ht7])

% Decide zoom level again if no flags were supplied, hide some labels this time
%if (nargin<1)
    if (cvg<scl(1))  % More zoom-in, probably half of a state
        tmp = findobj(VIS.fh,'String','County Seat');
        if ~isempty(tmp), figure(VIS.fh), eval(get(tmp(1),'Callback')); end
    end
    if (cvg<scl(2))  % Even more zoom-in
        tmp = findobj(VIS.fh,'String','Pop>10k');
        if ~isempty(tmp), figure(VIS.fh), eval(get(tmp(1),'Callback')); end
        tmp = findobj(VIS.fh,'String','Pop>1000');
        if ~isempty(tmp), figure(VIS.fh), eval(get(tmp(1),'Callback')); end
   end
%end


% If user does not want VIS variable, clear it
if nargout<1, clear VIS; end
return



function [x pos] = ismember(test_ele, ele_set)
pos = [];
if iscell(test_ele) && iscell(ele_set)
    x = false(size(test_ele));
    pos = zeros(size(test_ele));
    for jdx = 1:length(test_ele)
        for idx = 1:length(ele_set)
            if strcmp(test_ele{jdx},ele_set{idx}), x(jdx) = 1; pos(jdx) = idx; end
        end
    end
elseif iscell(test_ele) && ~iscell(ele_set)
    x = false;
    for idx = 1:length(test_ele)
        if strcmp(test_ele{idx},ele_set), x(idx) = 1; end
    end
elseif ~iscell(test_ele) && iscell(ele_set)
    x = false;
    for idx = 1:length(ele_set)
        if strcmp(test_ele,ele_set{idx}), x = 1; pos = idx; end
    end
end
return


function [yn] = isoverlap(L,R)
if numel(L)~=4
    fprintf('The first input must be a 4-element vector only.\n');
    yn = -1;
    return
end
if size(R,1)~=4,
    if numel(R)==4,
        R = R(:);
    else
        fprintf('Sorry I can''t compute.\n');
        yn = -1;
        return
    end
end
yn = ( (L(1)>R(1,:))&(L(1)<R(2,:))&(L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(1)>R(1,:))&(L(1)<R(2,:))&(L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&(L(3)>R(3,:))&(L(3)<R(4,:)) )|...
     ( (L(2)>R(1,:))&(L(2)<R(2,:))&(L(4)>R(3,:))&(L(4)<R(4,:)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&(R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(1,:)>L(1))&(R(1,:)<L(2))&(R(4,:)>L(3))&(R(4,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&(R(3,:)>L(3))&(R(3,:)<L(4)) )|...
     ( (R(2,:)>L(1))&(R(2,:)<L(2))&(R(4,:)>L(3))&(R(4,:)<L(4)) );
return
