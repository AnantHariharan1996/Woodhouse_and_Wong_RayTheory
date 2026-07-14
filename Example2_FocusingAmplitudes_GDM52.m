% Demonstrate the calculation of amplitudes using the 
% linearized approach from Woodhouse and Wong, 1986;
% This version just preloads phase velocity map.

clear; clc; close all; 
STANFO = readtable('STATIONS_TA');
slat = STANFO.Var3;
slon = STANFO.Var4;

% In this example, we load the station locations and event locations for
% one earthquake at one period. Then, we feed this information into the
% amplitude prediction code.
Period=100;

% Let's asume the event with CMT ID 201501100205A
elat = -5.72000;
elon =68.34000;
currevid = num2str(201501100205);
 % call the function which calculates amplitudes and give it the station and event locations
[Amplist] = Calculate_FocusingAmps_WW86_AllinOne_PreloadModels(elon,elat,slon,slat,Period,'GDM52');
[lattrck2,lonttrck2] = track2(elat,elon,mean(slat),mean(slon));

figure(1)
scatter(slon,slat,50,Amplist,'filled','MarkerEdgeColor','k','LineWidth',2)
hold on
scatter(elat,elon,100,'pentagram','filled')
plot(lonttrck2,lattrck2,'linewidth',2,'linestyle','--','color','k')
colormap(flipud(polarmap))
barbar=colorbar;
ylabel(barbar,'Focusing-Predicted Amplitudes')
caxis([0.4 4.6])
xlim([-130 -60])
ylim([20 50])
title([ currevid ': Predictions at ' num2str(Period) 's, assuming GDM52 Maps' ])
box on;
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')

