function [predicted_ttimelist] = Calculate_TravelTime_AllinOne_ArbitraryModel(lon_src,lat_src,lon_stas,lat_stas,lon_c,lat_c,c)
% Calculate ray-theoretic traveltimes corresponding to a given model.
% You need to provide this function with 7 things:
% lon_src: The longitude of the Earthquake source
% lat_src: The latitude of the Earthquake source
% lon_stas: The longitudes of the stations, or coordinates at which you
% want to predict the amplitudes
% lat_stas:The latitudes of the stations, or coordinates at which you
% want to predict the amplitudes
% lon_c: Longitudes at which the phase velocity map is defined.
% lat_c: Latitudes at whcih the phase velocity map is defined.
% c: Phase velocity map (km/s)
% Anant Hariharan, 2026
distlist = distance(lat_src,lon_src,lat_stas,lon_stas);
predicted_ttimelist = zeros(size(distlist));
for ijk = 1:length(lon_stas)
100*ijk/length(lon_stas)
currlon = lon_stas(ijk);
currlat = lat_stas(ijk);
[lattrk,lonttrk] = track2(currlat,currlon,lat_src,lon_src,[],'degrees',400);
ctrk = griddata(lon_c,lat_c,c,lonttrk,lattrk);
meanc = mean(ctrk);
currdist = deg2km(distlist(ijk));
predicted_ttimelist(ijk) = currdist./meanc;

end



end