% OCES5300-Matlab script for BGC-Argo file
% Save as a text file (.txt) to allow its submission on Canvas
% 20810279_LIN-Yuxuan

%% Import data file

clear all; close all; clc;

% targeted varibles:
% 'TEMP_ADJUSTED','PSAL_ADJUSTED','DOX2_ADJUSTED','CPHL_ADJUSTED','NTAW_ADJUSTED'
Temp = ncread('GL_PR_PF_5906305.nc','TEMP_ADJUSTED');
Temp = Temp(:,1);
DO = ncread('GL_PR_PF_5906305.nc','DOX2_ADJUSTED');
DO = DO(:,1);
NO3 = ncread('GL_PR_PF_5906305.nc','NTAW_ADJUSTED');
NO3 = NO3(:,2);
Chl = ncread('GL_PR_PF_5906305.nc','CPHL_ADJUSTED');
Chl = Chl(:,1);
Sal = ncread('GL_PR_PF_5906305.nc','PSAL_ADJUSTED');
Sal = Sal(:,1);
Depth = ncread('GL_PR_PF_5906305.nc','PRES_ADJUSTED');
Depth = Depth(:,1);
NO3_Depth = ncread('GL_PR_PF_5906305.nc','PRES_ADJUSTED');
NO3_Depth = NO3_Depth(:,2);

%% Plot information

%---------------------
% Plot 1: Profiles
%---------------------
figure(1);

% Temperature
subplot(1,5,1);
plot(Temp, Depth,'color','k','linewidth',2);
xlabel('Temperature [°C]', 'fontsize', 14, 'fontweight','bold');
ylabel('Depth (m)', 'fontsize', 14, 'fontweight','bold');

% Salinity
subplot(1,5,2);
plot(Sal, Depth,'color','k','linewidth',2);
xlabel('Salinity [‰]', 'fontsize', 14, 'fontweight','bold');

% Chlorophyll-a
subplot(1,5,3);
plot(Chl, Depth,'color','k','linewidth',2);
xlabel('Chlorophyll-a [mg m^{-3}]', 'fontsize', 14, 'fontweight','bold');


% Nitrate
subplot(1,5,4);
% Negative value for the first 18 records. Manually modified by zero
NO3(NO3<0) = 0;
plot(NO3, NO3_Depth,'color','k','linewidth',2);
xlabel('Nitrate [µmol kg^{-1}]', 'fontsize', 14, 'fontweight','bold');

% Dissolved oxygen
subplot(1,5,5);
plot(DO, Depth,'color','k','linewidth',2);
xlabel('Dissolved Oxygen [µmol kg^{-1}]', 'fontsize', 14, 'fontweight','bold');

% Personalization
for i = 1:5
    subplot(1,5,i);
    set(gca,'linewidth', 1.5, 'fontsize',12, 'Ydir', 'reverse', 'XAxisLocation', 'Top', 'box', 'off','YLim', [0,2100]);
end

%---------------------
% Plot 2: T-S diagram
%---------------------
figure(2);

% Determine graph boundaries
smin=min(Sal)-0.01.*min(Sal);
smax=max(Sal)+0.01.*max(Sal);
tmin=min(Temp)-0.1*max(Temp);
tmax=max(Temp)+0.1*max(Temp);

% Create a matrix for density contour 
xdim=round((smax-smin)./0.1+1);
ydim=round((tmax-tmin)+1);
dens=zeros(ydim,xdim);

% Calculate one-atmosphere density of seawater
ti=((1:ydim)-1)*1+tmin;
si=((1:xdim)-1)*0.1+smin;
for j=1:ydim
    for i=1:xdim
        % Define function of density calculator (valid for practical salinity from 0 to 42 and temperature from –2 to 40°C)
        % Reference: Pilson 2013, p59
        % An alternative choice is to separately define a density calculator function. Here I merge it into a single script
        t = ti(j);s = si(i);
        
        % Density of the Standard Mean Ocean Water (SMOW) 
        dw = 999.842594 + 6.793952*10^(-2)*t - 9.09529*10^(-3)*t^2 + 1.001685*10^(-4)*t^3 - 1.120083*10^(-6)*t^4 + 6.536332*10^(-9)*t^5;

        % The One Atmosphere International Equation of State of Seawater
        dens(j,i) = dw + (8.24493*10^(-1) - 4.0899*10^(-3)*t + 7.6438*10^(-5)*t^2 - 8.2467*10^(-7)*t^3 + 5.3875*10^(-9)*t^4)*s ...
        + (-5.72466*10^(-3) + 1.0227*10^(-4)*t - 1.6546*10^(-6)*t^2)*s^(3/2) + 4.8314*10^(-4)*s^2;
    end
end

% Draw T-S diagram
dens=dens-1000;
[c,h]=contour(si,ti,dens,'k');
clabel(c,h,'LabelSpacing',1000);
xlabel('Salinity [‰]','FontWeight','bold','FontSize',12);
ylabel('Temperature [°C]','FontWeight','bold','FontSize',12);

% Add scatter
hold on;
plot(Sal,Temp,'color','k','linewidth',2.5);

% References: Libes 2009; Emery 2001; Sarmiento and Gruber 2008; Jeandel et al 2013
% Potential water masses at different depths
% (1) Upper waters (<500 m): 
    % Eastern South Pacific Central Water (ESPCW, 8-24 °C, 34.4-36.4 PSU)
% (2) Intermediate waters (500-1500 m):
    % Eastern South Pacific Intermediate Water (ESPIW, 10-12 °C, 34-34.4 PSU)
    % Antarctic Intermediate Water (AAIW, 2-10 °C, 33.8-34.5 PSU, 26.8-27.4 sigma);
% (3) Deep and abyssal waters (>1500 m):
    % Upper Circumpolar Deep Waters (UCDW, 1.5-1.9 °C, 34.64-34.69 PSU)

hold on;
for d=[100 200 500 1000 1500 2000]
    d_r = Depth(abs(Depth - d) == min(abs(Depth - d)));
    scatter(Sal(Depth==d_r), Temp(Depth==d_r),70,'k','|');
    if d==100 || d==200
        x = 0.04; y = -0.2;
    elseif d==500
        x = 0.04; y = -0.1;
    elseif d==1000 || d==1500 || d==2000
        x = 0.03; y = 0.5;
    end
    text(Sal(Depth==d_r)+x, Temp(Depth==d_r)+y, num2str(d),'FontSize',12);
end

% Eastern South Pacific Central Water (ESPCW, 8-24 °C, 34.4-36.4 PSU)
ecc = axes2ecc(1.8,0.2);
[elat,elon] = ellipse1(34.7,19,[1.8 ecc],85);
plot(elat,elon,'linewidth',2,'color',[0 0.4470 0.7410]);
text(34.95, 18.6, 'ESPCW','FontWeight','bold','FontSize',17,'color',[0 0.4470 0.7410]);

% Eastern South Pacific Intermediate Water (ESPIW, 10-12 °C, 34-34.4 PSU)
ecc = axes2ecc(1,0.1);
[elat,elon] = ellipse1(34.2,11,[1 ecc],90);
plot(elat,elon,'linewidth',2,'color',[0.8500 0.3250 0.0980]);
text(34.35, 11, 'ESPIW','FontWeight','bold','FontSize',17,'color',[0.8500 0.3250 0.0980]);

% Antarctic Intermediate Water (AAIW, 2-10 °C, 33.8-34.5 PSU, 26.8-27.4 sigma);
ecc = axes2ecc(2,0.2);
[elat,elon] = ellipse1(34.2,5,[2 ecc],90);
plot(elat,elon,'linewidth',2,'color',[0.9290 0.6940 0.1250]);
text(34.45, 6, 'AAIW','FontWeight','bold','FontSize',17,'color',[0.9290 0.6940 0.1250]);

% Upper Circumpolar Deep Waters (UCDW, 1.5-1.9 °C, 34.64-34.69 PSU)
ecc = axes2ecc(0.8,0.07);
[elat,elon] = ellipse1(34.65,2,[0.8 ecc],90);
plot(elat,elon,'linewidth',2,'color',[0.4940 0.1840 0.5560]);
text(34.75, 2, 'UCDW','FontWeight','bold','FontSize',17,'color',[0.4940 0.1840 0.5560]);

%% Display summary

% === Q1: Location of float cast === %
% Can also read from attributes: Lat = ncreadatt('GL_PR_PF_5906305.nc','/','geospatial_lat_min');
% 'min' is equal to 'max'. Same as Longitude
% Same as 'Last observation'
Lat = ncread('GL_PR_PF_5906305.nc','LATITUDE');
Lat = Lat(1,:);
Lon = ncread('GL_PR_PF_5906305.nc','LONGITUDE');
Lon = Lon(1,:);

% === Q2: Time of float cast === %
% Can also read from attributes: time = ncreadatt('GL_PR_PF_5906305.nc','/','time_coverage_start');
% 'start' is equal to 'end' 
% Same as 'Last observation'
time = datestr(ncread('GL_PR_PF_5906305.nc','TIME')+712224);
time = time(1,1:20);

% === Q3: Depth of the surface mixed layer based on (static) temperature === %
% A temperature change threshold from the ocean subsurface of 0.2°C (i.e., Zref=10 m, deltaT=0.2°C) was chosen (de Boyer Montégut et al. 2004; Libes 2009, p70)
% Sensitivity: slightly deeper MLD of 49.94 m based on Monterey and Levitus (1997) (i.e., Zref=0 m, deltaT=0.5°C, details see: https://www.nodc.noaa.gov/OC5/WOA94/mix.html)
Zref = Depth(abs(Depth - 10) == min(abs(Depth - 10)));
MLD = Depth(find(Temp < Temp(Depth==Zref)-0.2,1)-1);

% === Q4: Average surface salinity of the surface mixed layer === %
Sali_ML = mean(Sal(Depth <= MLD));

% === Q5: Depth of the deep chlorophyll maximum and its chlorophyll concentration === %
DCM = Depth(Chl==max(Chl));
Chl_DCM = max(Chl);

% === Q6: Average DO in the surface mixed layer relative to the average DO below that to 200 m depth === %
DO_ML = mean(DO(Depth <= MLD));
DO_below = mean(DO(Depth > MLD & Depth <= 200));
R_DO = DO_ML/DO_below;

% === Q7: Depth and oxygen concentration at the deep oxygen minimum === %
d_OMZ = Depth(DO == min(DO));
OMZ = min(DO);

% === Q8: Nitrate at the depth of the oxygen minimum compared to the average nitrate in the mixed layer === %
% Nodata for nitrate at OMZ
% Negative nitrate at surface layer
d_OMZ_NO3 = NO3_Depth(abs(NO3_Depth - d_OMZ)==min(abs(NO3_Depth - d_OMZ)));
NO3_OMZ = NO3(NO3_Depth == d_OMZ_NO3);
NO3_ML = mean(NO3(NO3_Depth <= MLD));
delta_NO3 = NO3_OMZ - NO3_ML;

% === Q9: Number of conservative water masses undergoing conservative mixing === %
% See plot 2 for water mass identification
nwm = 4;

% QA summaries
disp(['___________________________']) 
disp([' ']) 
disp(['20810279_LIN_Yuxuan data summary / answers:'])
disp([' ']) 
disp(['Q1: The location of float cast (in degree) = Latitude ' num2str(Lat) '; Longitude ' num2str(Lon) '; in Southeastern Pacific Ocean'])
disp(['------']) 
disp(['Q2: Time of float cast = ' num2str(time)])
disp(['------']) 
disp(['Q3: Mixed layer depth = ' num2str(MLD) ' [m]'])
disp(['------']) 
disp(['Q4: Average salinity of the mixed layer = ' num2str(Sali_ML) ' [‰]'])
disp(['------']) 
disp(['Q5: Depth of the deep chlorophyll maximum = ' num2str(DCM) ' [m]; the DCM concentration = ' num2str(Chl_DCM) ' [mg m-3]'])
disp(['------']) 
disp(['Q6: Average DO in the surface mixed layer relative to the average DO below that to 200 m depth = ' num2str(R_DO)])
disp(['------']) 
disp(['Q7: Depth of the deep oxygen minimum = ' num2str(d_OMZ) ' [m] ; the minimum oxygen concentration = ' num2str(OMZ) ' [µmol kg-1]'])
disp(['------']) 
disp(['Q8: Nitrate at the depth of the oxygen minimum compared to the average nitrate in the mixed layer = ' num2str(delta_NO3)])
disp(['------']) 
disp(['Q9: Number of conservative water masses undergoing conservative mixing = ' num2str(nwm)])
disp([' #1 Eastern South Pacific Central Water (ESPCW, <500 m)'])
disp([' #2 Eastern South Pacific Intermediate Water (ESPIW, <500 m)'])
disp([' #3 Antarctic Intermediate Water (AAIW, 500-1500 m)'])
disp([' #4 Upper Circumpolar Deep Waters (UCDW, >1500 m)'])
disp(['Water masses are identified in Plot 2 in colored font and lines'])
disp([' ']) 
disp(['___________________________']) 

