% Demonstrate the calculation of arrival angles using the 
% linearized approach from Woodhouse and Wong, 1986;
% This is sometimes referred to as the 'path integral approximation'. 

clear; clc; close all; 
EVIDLIST = {'200903122323A'};

% In this example, we load the station locations and event locations for
% one earthquake at one period. Then, we feed this information into the
% arrival angle prediction code and it generates theoretical arrival
% angles, and then plots them.
Period=100;
for ijk =1:length(EVIDLIST)
   currevid = EVIDLIST{ijk};


    fname = ['GoodAmpdtpamp_50s_' currevid ];
    info = load(fname);
% load station locations and the location of the earthquake
    slat = info(:,3); % station latitude
    slon = info(:,4); % station longitude
    elat = info(1,1); %earthquake latitude
    elon = info(1,2); %earthquake longitude

    
    % call the function which calculates arrival angles.
    % give this function the station and event locations


    % Get the predictions for the GDM52 maps!
[AAlist] = Calculate_Arrival_Angle_WW86_AllinOne_PreloadModels(elon,elat,slon,slat,Period,'GDM52');

figure(1)
subplot(1,2,1)
scatter(slon,slat,50,AAlist,'filled','MarkerEdgeColor','k','LineWidth',2)
colormap(flipud(polarmap))
barbar=colorbar;
caxis([-20 20])
xlim([-125 -85])
ylim([25 50])
title([currevid ': Predictions at ' num2str(Period) 's, GDM52 Maps' ])
box on;
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')
ylabel(barbar,'Arrival Angle (deg)')
% Now, do the same thing, but for Zhitu's phase velocity map

[AAlist2] = Calculate_Arrival_Angle_WW86_AllinOne_PreloadModels(elon,elat,slon,slat,Period,'Ma2014');


subplot(1,2,2)
scatter(slon,slat,50,AAlist,'filled','MarkerEdgeColor','k','LineWidth',2)
colormap(flipud(polarmap))
barbar=colorbar;
caxis([-20 20])
xlim([-125 -85])
ylim([25 50])
title([currevid ': Predictions at ' num2str(Period) 's, Ma 2014 Maps' ])
box on;
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')
ylabel(barbar,'Arrival Angle (deg)')

set(gcf,'position',[110 388 1269 418])

%%%%



figure(100)
scatter(AAlist,AAlist2,50,'filled')
hold on
xlim([-10 10])
ylim([-10 10])
plot([-10 10],[-10 10],'linewidth',2,'color','k','linestyle','--')
box on
title('Comparison Between Predictions from two Models')
xlabel('GDM52 Predictions (deg)')
ylabel('Ma 2014 Predictions(deg)')
set(gca,'fontsize',18,'linewidth',2,'fontweight','bold')



end