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

data_file=[data_path deployment '_ocean_glider_velo.mat'];
load(data_file);

%% Split the ocean velocity profiles into E/W and N/S components
% u=real(O_ls) to get real (E/W) component
% v=imag(O_ls) to get imag (N/S) component
u=real(O_ls);
v=imag(O_ls);

%% Calculate shear from ocean velocity
% Find min and max depth range for the grid from bnew
d_min=min(bnew);
d_max=max(bnew);
dz   =bnew(2)-bnew(1); % Vertical resolution of the profiles

% Find the necessary number of depth bins
bin_num=((d_max-d_min)./dz)+1; % Add 1 to account for d_min bin

% Number of profile columns needed to fill with ocean velocity data
% Each profile starts with a depth less than the depth recorded at the end
% of the last profile.
prof_num=sum(diff(bnew)<0)+1;

% Make a grid to populate with shear values
% Create the grid
usgrid=nan([bin_num prof_num]);
vsgrid=nan([bin_num prof_num]);
s_maggrid=nan([bin_num prof_num]); % Grid for shear magnitude

% Create a reference array for bin depth index
grid_bin=(d_min:dz:d_max)';

% Loop through each profile and find the du/dz dv/dz and then populate a
% grid
profile_ind=[0; find(diff(bnew)<0)]+1; % Index of the start of each profile

for ii=1:length(profile_ind)-1
    if ~isnan(profile_ind(ii))
        % Grabs one profile of data
        prof_bins=bnew(profile_ind(ii):profile_ind(ii+1)-1)';
        prof_u   =u(profile_ind(ii):profile_ind(ii+1)-1)';
        prof_v   =v(profile_ind(ii):profile_ind(ii+1)-1)';
        
        % Calculate the shear in each direction
        prof_us=diff(prof_u)./diff(prof_bins);
        prof_vs=diff(prof_v)./diff(prof_bins);
        prof_mags=sqrt(prof_us.^2 + prof_vs.^2); % Shear magnitude
        
        % Finds the depths to at which to insert this data into the grid
        [prof_mindep, grid_start]=min(prof_bins);
        [prof_maxdep, grid_end]=max(prof_bins);
    
        % Populate the grid with shear calculations
        % Need to subtract one off the index because diff returns N-1 points
        usgrid(grid_start:grid_end-1,ii)=prof_us;
        vsgrid(grid_start:grid_end-1,ii)=prof_vs;
        s_maggrid(grid_start:grid_end-1,ii)=prof_mags;
    else
    end
end
save([data_path deployment '_ocean_shear.mat'], 'usgrid', 'vsgrid', 's_maggrid', 'grid_bin');

%% Plotting gridded shear in each direction

% Get some colorbar limits. Define limits the same for u and v
tolerance=3; % standard deviation tolerance for the colorbar limits
shear=[usgrid(:);vsgrid(:)];
mshear=nanmean(shear);
stdshear=nanstd(shear);
clim=[mshear-tolerance.*stdshear mshear+tolerance.*stdshear];

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