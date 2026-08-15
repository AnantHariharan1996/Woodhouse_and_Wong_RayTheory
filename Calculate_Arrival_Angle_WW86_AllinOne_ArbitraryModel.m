function [AAlist] = Calculate_Arrival_Angle_WW86_AllinOne_ArbitraryModel(lon_src,lat_src,lon_stas,lat_stas,lon_c,lat_c,c)
%%% Full Implementation of woodhose and wong's linearized approach to
%%% predict arrival angles for any source-receiver path. 
% You need to provide this function with 7 things:
% lon_src: The longitude of the Earthquake source
% lat_src: The latitude of the Earthquake source
% lon_stas: The longitudes of the stations, or coordinates at which you
% want to predict the amplitudes
% lat_stas:The latitudes of the stations, or coordinates at which you
% want to predict the amplitudes
% lon_c: Longitudes at which the phase velocity map is defined.
% lat_c: Latitudes at whcih the phase velocity map is defined.
% c: Phase velocity map (km/s or m/s; since this is eventually non-dimensionalized, the units don't matter). 
% Anant Hariharan, 2025

%preallocate output
AAlist = NaN*ones(size(lon_stas));


lon = lon_c;
lat = lat_c;

[lat_TEMP, lon_TEMP] = fibonacci_sphere(64800);
c_TEMP = griddata(lon,lat,c,lon_TEMP,lat_TEMP);
c_ref = nanmean(c_TEMP(:));

% discretization along path and perpendicular to path
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




function [lat, lon] = fibonacci_sphere(samples)
    % Generate points on a sphere using the Fibonacci method
    % Output in degrees: latitude (-90 to 90) and longitude (-180 to 180)

    if nargin < 1
        samples = 1000; % default
    end

    phi = pi * (sqrt(5) - 1); % golden angle in radians

    lat = zeros(samples, 1);
    lon = zeros(samples, 1);

    for i = 0:samples-1
        y = 1 - (i / (samples - 1)) * 2; % y from 1 to -1
        radius = sqrt(1 - y^2);          % radius at y

        theta = phi * i;

        x = cos(theta) * radius;
        z = sin(theta) * radius;

        % Convert to latitude/longitude
        lat(i+1) = asind(y);                   % latitude in degrees
        lon(i+1) = atan2d(z, x);               % longitude in degrees
    end
end