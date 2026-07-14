function [AAlist] = Calculate_Arrival_Angle_WW86_AllinOne_PreloadModels(lon_src,lat_src,lon_stas,lat_stas,Period,Model)
%%% Full Implementation of woodhose and wong's linearized approach to
%%% predict arrival angles for any source-receiver path. 
%%% Right now, this is just set up to work with GDM52 or zhitu's phase velocity maps, 
%%% so you need to have those  files in your path! 
%
% You need to provide this function with 5 things:
% lon_src: The longitude of the Earthquake source
% lat_src: The latitude of the Earthquake source
% lon_stas: The longitudes of the stations, or coordinates at which you
% want to predict the arrival angles
% lat_stas:The latitudes of the stations, or coordinates at which you
% want to predict the arrival angles
% Period: The period you are interested in. This tells the code which phase
% velocity map to load. 
% June 29: This version of the code
% can choose between different models; the GDM52 maps
% of Ekstrom 2011 and Ma 2014. 
% Specify 'GDM52' or 'Ma2014'
% Anant Hariharan, 2025

if strcmp(Model,'GDM52')
% load the phase velocity maps 
PeriodStr = num2str(Period);
if length(PeriodStr) < 3
PeriodStr=['0' PeriodStr];
end
GDMFname = ['R' PeriodStr '_1_GDM52.pix'];
if exist(GDMFname) ~= 2
error('Looks like the Phase Velocity map from Ekstrom et al., 2011 does not exist at this frequency. You need to use a different period.')
end


 [ lon,lat,phvel,phvel_ref ] = Read_GDM52_Phvels( GDMFname );
lon=lon;
lat=lat;
c=phvel;
c_ref=phvel_ref;



elseif strcmp(Model,'Ma2014')
% load the phase velocity maps 
freq_mhz = 1000/Period;
freq_mhz_str = num2str(freq_mhz);
if freq_mhz < 10
freq_mhz_str = ['0' num2str(freq_mhz)];
end
MaFname = ['map_rayleigh_' freq_mhz_str '.txt'];
if exist(MaFname) ~= 2
error('Looks like the Phase Velocity map from Ma et al., 2014 does not exist at this frequency. You need to use a different period.')
end
info = load(MaFname);
lon = info(:,1);
lat =info(:,2);
c = info(:,3);
c_ref  =mean(c);
end




% 
N=500;
integral_spacing = 0.1;

% First, get the perturbations wrt reference phvel.
dcc = (c-c_ref)./c_ref;

% Loop over stations and get the Arrival angle for every station
for stanum = 1:length(lon_stas)

disp([ 'Completed: ' num2str(100*stanum/length(lon_stas) ) '% of total stations for this event'  ])
lon_sta = lon_stas(stanum);
lat_sta =lat_stas(stanum);

% rottate to equator
 [lon_rot,lat_rot] = greatcircle_fast(lon,lat,lon_src,lat_src,lon_sta,lat_sta);
[staxrot,stayrot] = greatcircle_fast(lon_sta,lat_sta,lon_src,lat_src,lon_sta,lat_sta);

% get points along the path. 
[~,lon_gc] = track2(0,0,stayrot,staxrot,[],'degrees',N);

% get points just above the path
upperlats = integral_spacing*ones(size(lon_gc));
% get points just below the path.  
% Remember, this is the equator...
lowerlats =-1*integral_spacing*ones(size(lon_gc));

% As above so below; get the phase velocities at these coordinates
Values_Upper = griddata(lon_rot(:),lat_rot(:),dcc(:),lon_gc,upperlats);
Values_Lower = griddata(lon_rot(:),lat_rot(:),dcc(:),lon_gc,lowerlats);
difference_wrt_colatitude = Values_Lower-Values_Upper;
dtheta = deg2rad(integral_spacing);
% finite difference derivative (centered?)
deriv_wrt_colat = difference_wrt_colatitude/(2*dtheta);
% now perform integration

% get to radians
longc_rad = deg2rad(lon_gc);

% Thanks, woodhouse
integrand = sin(longc_rad).*deriv_wrt_colat;
Delta  = distance(lat_src,lon_src,lat_sta,lon_sta);
zeta = -1/sin(   deg2rad(Delta)    );
zeta=zeta*trapz(longc_rad,integrand);
% get back to degrees. 
AA = rad2deg(zeta);


AAlist(stanum) = AA;
end

end



function [lon_rot,lat_rot] = greatcircle_fast(lon,lat,lon_src,lat_src,lon_sta,lat_sta)
% rotate everything to the great-circle path
% radians
lon = deg2rad(lon); lat = deg2rad(lat);
lon_src = deg2rad(lon_src); lat_src = deg2rad(lat_src);
lon_sta = deg2rad(lon_sta); lat_sta = deg2rad(lat_sta);

% source and receiver unit vectors
s = [cos(lat_src)*cos(lon_src);
     cos(lat_src)*sin(lon_src);
     sin(lat_src)];

r = [cos(lat_sta)*cos(lon_sta);
     cos(lat_sta)*sin(lon_sta);
     sin(lat_sta)];

% basis vectors
ex = s;
ez = cross(s,r); ez = ez/norm(ez);
ey = cross(ez,ex);

% convert all points to Cartesian
x = cos(lat).*cos(lon);
y = cos(lat).*sin(lon);
z = sin(lat);

% projections (fast!)
xp = ex(1)*x + ex(2)*y + ex(3)*z;
yp = ey(1)*x + ey(2)*y + ey(3)*z;
zp = ez(1)*x + ez(2)*y + ez(3)*z;

% back to spherical
lat_rot = rad2deg(atan2(zp,sqrt(xp.^2+yp.^2)));
lon_rot = rad2deg(atan2(yp,xp));

end

function [ lon,lat,phvel,phvel_ref ] = Read_GDM52_Phvels( fname )
% Reads GDM52 Phase velocity files.
% These files are downloaded from the GDM52 website
% https://www.ldeo.columbia.edu/~ekstrom/Projects/SWP/GDM52.html
% These files  must be in the same folder as this function. 
% Read reference phase velocity
fid=fopen(fname);
linenum=3;
C=textscan(fid,'%s %f',1,'delimiter','\n','headerlines',linenum-1);
GrpvelString=C{1};
GrpvelString=GrpvelString{1};
GrpvelString=GrpvelString(end-7:end);
fclose(fid);
referencegrpvel=str2num(GrpvelString);
% Read information
fid=fopen(fname);
data=textscan(fid,'%f %f %f %f','HeaderLines',6);
fclose(fid);

lat=data{1};
lon=data{2};
grpveldeviation=data{4};

phvel=referencegrpvel+referencegrpvel*grpveldeviation/100;

phvel_ref=referencegrpvel;
end

