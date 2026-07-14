% Demonstrate the focusing of amplitudes that results from a low-velocity
% anomaly. First, generate a synthetic phase velocity map with a
% low-velocity anomaly. Then use the expressions in Woodhouse and Wong 1986
% to predict the amplitudes that result, 
% assuming an earthquake located at -20,20
% Schematically, this should look pretty similar to the diagram 
% in Feng and Ritzwoller, GJI

clear; clc; close all; 
load coastlines

%%% Step 1: Specify the input phase velocity
% Very crude way to generate synthetic phase velocity map with
% a single focusing blobby low velocity anomaly
LonList = [-179.5:0.75:179.5];
LatList = [-89.5:0.75:89.5];
[lon_c,lat_c] = meshgrid(LonList,LatList);
CapRadius = 3; % degrees; radius of spherical cap of low velocity anomaly
Centroid_Lon = 30;
Centroid_Lat = 20;
dist2pts = distance(Centroid_Lat,Centroid_Lon,lat_c,lon_c);
% Make it kinda gaussian
gaussnfn = -0.1*exp(-0.5 * (( - dist2pts) /2.5).^2);
c = 4+0.75*(gaussnfn);
lon_c=lon_c(:); lat_c=lat_c(:); c=c(:);

%%% Step 2: specify station locations
Sta_X = [20:2.5:47.5];
Sta_Y = [15:2.5:50];
[lon_stas,lat_stas]=meshgrid(Sta_X,Sta_Y);
lon_stas=lon_stas(:);
lat_stas=lat_stas(:);

%%% Step 3: Specify event location
lat_src = -20;
lon_src =20;
%%%%%%%%%
%%% You should not need to change anything below here

tic
 % call the function which calculates amplitudes and give it the station and event locations
[Amplist] = Calculate_FocusingAmps_WW86_AllinOne_ArbitraryModel(lon_src,lat_src,lon_stas,lat_stas,lon_c,lat_c,c);
toc

[lattrck2,lonttrck2] = track2(lat_src,lon_src,Centroid_Lat,Centroid_Lon);



figure()
subplot(1,2,1)
scatter(lon_c,lat_c,25,c,'filled')
hold on
scatter(lon_stas,lat_stas,10,[1 0 1],'^','filled')
scatter(lon_src,lat_src,200,[1 1 1],'pentagram','filled')
hold on
plot(coastlon,coastlat,'linewidth',2,'color','k')
barbar=colorbar;
ylabel(barbar,'Phase Velocity (km/s)')
title('Synthetic Phase Velocity Map with Low Velocity Anomaly')
xlim([-30 80])
ylim([-60 60])
box on;
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')


subplot(1,2,2)
scatter(lon_stas,lat_stas,35,Amplist,'filled','MarkerEdgeColor','k','LineWidth',0.5)
hold on
scatter(lon_src,lat_src,200,'pentagram','filled')
plot(lonttrck2,lattrck2,'linewidth',2,'linestyle','--','color','k')
colormap(flipud(polarmap))
barbar=colorbar;
plot(coastlon,coastlat,'linewidth',2,'color',[0.5 0.5 0.5])
ylabel(barbar,'Focusing-Predicted Amplitudes')
xlim([-30 80])
ylim([-60 60])
caxis([0.8 1.2])
colormap(flipud(jet))
title(['Amplitudes due to Focusing from Low-Velocity Anomaly'])
box on;
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')
set(gcf,'position',[-17 241 1451 579])

saveas(gcf,'Example1.png')
