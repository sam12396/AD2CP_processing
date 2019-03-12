%% Calculating shear profiles

%% Notes
% This program takes velocity profile data that was output by 
% AD2CP_ls_inversion.m and calculates the shear profile of each of the velocity
% profile. 

% Sam Coakley 2/29/2019

%% Starters
clear all; close all
deployment='ru31_467_2016_01_13';
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment '/'];
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' deployment '/'];

% Load the velcoity data to calculate shear
data_file=[data_path deployment '_ocean_velo_grid.mat'];
load(data_file);

%% Calculate shear from ocean velocity

% Get dimensions of grid
bin_num=size(ugrid,1); 
prof_num=size(ugrid,2);
dz =grid_bin(2)-grid_bin(1);
% Make a grid to populate with shear values
% Create the grid
usgrid=nan([bin_num prof_num]);
vsgrid=nan([bin_num prof_num]);
s_maggrid=nan([bin_num prof_num]); % Grid for shear magnitude

% Loop through each profile and find the du/dz dv/dz and then populate a
% grid
for ii=1:prof_num
    prof_us=diff(ugrid(:,ii))./diff(grid_bin);
    prof_vs=diff(vgrid(:,ii))./diff(grid_bin);
    prof_mags=sqrt(prof_us.^2 + prof_vs.^2); % Shear magnitude
    
    usgrid(1:end-1,ii)=prof_us;
    vsgrid(1:end-1,ii)=prof_vs;
    s_maggrid(1:end-1,ii)=prof_mags;

end
save([data_path deployment '_ocean_shear.mat'], 'usgrid', 'vsgrid', 's_maggrid', 'grid_bin');

%% Plotting gridded shear in each direction

% Get some colorbar limits. Define limits the same for u and v
tolerance=3; % standard deviation tolerance for the colorbar limits
shear=[usgrid(:);vsgrid(:)];
mshear=nanmean(shear);
stdshear=nanstd(shear);
clim=[-(mshear+tolerance.*stdshear) mshear+tolerance.*stdshear]; % center color bar at 0

% Set plotting constants
cmap=cmocean('delta');

% Plotting
figure(1)
pcolorjw(1:prof_num,-grid_bin,usgrid)
ax=gca;
hold on
cb=colorbar;
cb_label=get(cb,'Label');
set(cb_label,'String','East (+)/ West(-) Vertical Shear [1/s]');
ax.CLim=clim;
xlabel('Segment number')
ylabel('Depth [m]')
colormap(cmap)
title([deployment ': East (+)/ West(-) Vertical Shear [1/s] Vertical resolution:' num2str(dz) 'm'],'Interpreter','none')
print([fig_path 'u_ocean_shear_cs'],'-dpng')
close all
clear ax cb cb_label

figure(2)
pcolorjw(1:prof_num,-grid_bin,vsgrid)
ax=gca;
hold on
cb=colorbar;
cb_label=get(cb,'Label');
set(cb_label,'String','North (+)/ South(-) Vertical Shear [1/s]');
ax.CLim=clim;
xlabel('Segment number')
ylabel('Depth [m]')
colormap(cmap)
title([deployment ': North (+)/ South(-) Vertical Shear [1/s] Vertical resolution:' num2str(dz) 'm'],'Interpreter','none')
print([fig_path 'v_ocean_shear_cs'],'-dpng')
close all
clear ax cb cb_label

%% Plot gridded shear magnitude

% Get some colorbar limits. Define limits the same for u and v
tolerance=3; % standard deviation tolerance for the colorbar limits
shear=s_maggrid(:);
mshear=nanmean(shear);
stdshear=nanstd(shear);
clim=[0 mshear+tolerance.*stdshear];

% Set number of contours
contour_num=5;
contour_interval=clim(2)/contour_num;

% Define colormap
cmap=cmocean('speed');

figure(1)
contourf(1:prof_num,-grid_bin,s_maggrid,clim(1):contour_interval:clim(2))
ax=gca;
hold on
cb=colorbar;
cb_label=get(cb,'Label');
set(cb_label,'String','Magntiude of Vertical Shear [1/s]');
ax.CLim=clim;
xlabel('Segment number')
ylabel('Depth [m]')
colormap(cmap)
title([deployment ': Magntiude of Vertical Shear [1/s] Vertical resolution:' num2str(dz) 'm'],'Interpreter','none')
print([fig_path 'mag_ocean_shear_cs'],'-dpng')
close all
clear ax cb cb_label