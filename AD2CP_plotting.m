%% Plotting of t,s,d and velocity data calculated from AD2CP
% This script will make two plots. (1) Comparing t,s,d over the duration of
% the glider's mission. (2) Comparing Shear Magnitude, N2, and Richardson
% number over the same duration. The shear magnitude and Richardson number
% are bin x segment, N2 is bin x profile, and tsd is averaged into the depth
% bins. I have already made most of these plots individually in different 
% programs but putting it all together here will help with formatting.

% Sam Coakley
% 3/6/19

%% Starters
clear all;close all;
% Define the deployment of interest
deployment='ru31_467_2016_01_13';
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' deployment '/'];
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment '/'];

shear_data=load([data_path deployment '_ocean_shear.mat']);
g_data=load([data_path deployment '_dgroup.mat']);
ri_data=load([data_path deployment '_rich_num_grid.mat']);

%% Set some plotting parameters
std_threshold=3;
siz_title = 20;
siz_text = 16;
siz_marker= 8;

%% Load glider data
% Put the data to a grid with the same depth interval as the AD2CP data
[g_time, g_bin, temp]=g_data.dgroup.toGrid2d('drv_sea_water_temperature','depthbin',4);
[~, ~, sal]          =g_data.dgroup.toGrid2d('drv_sea_water_salinity','depthbin',4);
[~, ~, dens]         =g_data.dgroup.toGrid2d('drv_sea_water_density','depthbin',4);

%% Plotting glider t,s,d
% Make colorbar limits for each variable
tm=nanmean(temp(:));
ts=nanstd(temp(:));
tclim=[(tm-std_threshold*ts) (tm+std_threshold*ts)];
sm=nanmean(sal(:));
ss=nanstd(sal(:));
sclim=[(sm-std_threshold*ss) (sm+std_threshold*ss)];
dm=nanmean(dens(:));
ds=nanstd(dens(:));
dclim=[(dm-std_threshold*ds) (dm+std_threshold*ds)];

% Define colormaps for each variable
tcmap=cmocean('thermal');
scmap=cmocean('haline');
dcmap=cmocean('dense');

figure(1)
hold on

% Plot temperature
subplot(3,1,1)
pcolorjw(g_time,-g_bin,temp)
shading interp
tax=gca;
tcb=colorbar;
tax.CLim=tclim;
tax.Colormap=tcmap;
datetick('x',6)
ylabel('Depth [m]')
tcb_label=get(tcb,'Label');
set(tcb_label,'String','°C');
title([deployment ' Temperature, Salinity, and Density'],'Interpreter','none')

% Plot salinity
subplot(3,1,2)
pcolorjw(g_time,-g_bin,sal)
shading interp
sax=gca;
scb=colorbar;
sax.CLim=sclim;
sax.Colormap=scmap;
datetick('x',6)
ylabel('Depth [m]')
scb_label=get(scb,'Label');
set(scb_label,'String','PSU');

% Plot salinity
subplot(3,1,3)
pcolorjw(g_time,-g_bin,dens)
shading interp
dax=gca;
dcb=colorbar;
dax.CLim=dclim;
dax.Colormap=dcmap;
datetick('x',6)
ylabel('Depth [m]')
dcb_label=get(dcb,'Label');
set(dcb_label,'String','kg m^-^3');


%% Grab the shear and richardson data from the structures
s_maggrid=shear_data.s_maggrid;
grid_bin=shear_data.grid_bin;
N2_grid=ri_data.N2_grid;
N2_time=ri_data.profile_time;
ri_grid=ri_data.ri_grid;
segs=ri_data.segment_start;

%% Plotting shear, N2, and rich num
% Shear magnitude and rich num are on the same grid but N2 is on a bin x
% profile grid.
% Make colorbar limits for each variable
n2m=nanmean(N2_grid(:));
n2s=nanstd(N2_grid(:));
n2clim=[0 (n2m+std_threshold*n2s)];
smm=nanmean(s_maggrid(:));
sms=nanstd(s_maggrid(:));
smclim=[0 (smm+std_threshold*sms)];
% Rm=nanmean(ri_grid(:));
% Rs=nanstd(ri_grid(:));
% Rclim=[0 (Rm+std_threshold*Rs)];

% Define colormaps for each variable
n2cmap=cmocean('speed');
smcmap=cmocean('speed');
Rcmap=cmocean('speed');


figure(2)
subplot(3,1,1)
pcolorjw(N2_time,-grid_bin,N2_grid)
shading interp
n2ax=gca;
n2cb=colorbar;
n2ax.CLim=n2clim;
n2ax.Colormap=n2cmap;
datetick('x',6)
ylabel('Depth [m]')
n2cb_label=get(n2cb,'Label');
set(n2cb_label,'String','s^-^1');
title([deployment ' N2'],'Interpreter','none')

subplot(3,1,2)
pcolorjw(segs,-grid_bin,s_maggrid)
shading interp
smax=gca;
smcb=colorbar;
smax.CLim=smclim;
smax.Colormap=smcmap;
datetick('x',6)
ylabel('Depth [m]')
smcb_label=get(smcb,'Label');
set(smcb_label,'String','s^-^1');
title([deployment ' dU/dz'],'Interpreter','none')


subplot(3,1,3)
pcolorjw(segs,-grid_bin,ri_grid)
%contourf(segs,-grid_bin,ri_grid,[0 1 2 3])
hold on
shading interp
contour(segs,-grid_bin,ri_grid,[0.25 0.25],'color','m','LineWidth',1.5)
Rax=gca;
Rcb=colorbar;
Rax.CLim=[0 5];
Rax.Colormap=Rcmap;
datetick('x',6)
ylabel('Depth [m]')
Rcb_label=get(Rcb,'Label');
set(Rcb_label,'String','[no dim]');
title([deployment ' Ri'],'Interpreter','none')



