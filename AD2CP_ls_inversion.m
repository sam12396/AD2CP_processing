%% Splitting the velocity into ocean and glider components
% This program will use the function inversion_leastSquare_sparse_2019.m to
% separate the ocean velocity and glider velocity components from the AD2CP
% measured velocities. This will return one ocean velocity profile and one
% glider velocity profile for each segment of the glider's mission. 
% This is then used to calculate the vertical shear in each profile in each
% direction, (du/dz and dv/dz). 

% This is step 3 in our analysis of AD2CP data. After the inversion is
% complete, the data is ready for any type of analysis having been
% seperated into ocean and glider components at depths of the corresponding
% bins. Our goal is to calculate shear profiles from the data so we proceed
% to that step after the inversion.

% Sam Coakley
% 2/27/2019

%% Starters
clear all; close all
deployment='ru31_467_2016_01_13';
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment '/'];
g_path=[data_path deployment '_dgroup.mat'];
ad2cp_path=[data_path deployment '_AD2CP.nc'];
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' deployment '/'];

% Switch to run the inversion or not. This is to save time if you want to
% plot but do not need to re-run the data analysis.
if isfile([data_path deployment '_ocean_glider_velo.mat'])
    run_inv=0;
else
    run_inv=1;
end

%% Load in glider data
load(g_path);
% Get the start and end times for each segment
seg_start=dgroup.startDatenums;
seg_end  =dgroup.endDatenums;

% Get the dead reckoned glider velocity for each segment
% There are more values of velocity than there are segments so average them
% together? This is not publication ready data analysis
data=toArray(dgroup,'sensors',{'m_water_vx','m_water_vy'});
mean_water_vx=nan(length(seg_start),1);
mean_water_vy=nan(length(seg_start),1);
for jj=1:length(seg_start)
    segs=seg_start(jj);
    sege=seg_end(jj);
    seg_data= data(:,1)>= segs & data(:,1)<sege ;
    mean_water_vx(jj)=nanmean(data(seg_data,3));
    mean_water_vy(jj)=nanmean(data(seg_data,4));
end

%% Load in AD2CP data for this segment
% Get number of bins in AD2CP profiles
bin_num=ncread(ad2cp_path,'Average_Range');
bin_num=size(bin_num,2);

% Load AD2CP time to use it as an index to pull out just one segment
a_time=ncread(ad2cp_path,'Average_MatlabTimeStamp');

if run_inv==1 % Only run this loop if the data has not already been processed
    O_ls=[];
    G_ls=[];
    bnew=[];
    C=[];
    for jj=1:length(seg_start) % Loop through each segment
    % Find the AD2CP data inside a segment
    seg_time_ind= a_time>=seg_start(jj) & a_time<seg_end(jj);

    % Find the index window of data that is contained in segment
    seg_time_start=find(seg_time_ind,1,'first');
    seg_time_end=find(seg_time_ind,1,'last');
    seg_count=length(seg_time_start:seg_time_end);
    
    % E/W velocity
    u=ncread(ad2cp_path,'Average_EVelocity',[seg_time_start 1],[seg_count bin_num]); 
    % N/S velocity
    v=ncread(ad2cp_path,'Average_NVelocity',[seg_time_start 1],[seg_count bin_num]); 
    % Depth of measurements
    z=ncread(ad2cp_path,'Average_Depth',[seg_time_start 1],[seg_count bin_num]); 

    if jj==1 
        % How do we want to set dz to? 
        dz=ceil(z(1,2)-z(1,1)).*2; %Desired resolution of final profile
    end
    
    [tO_ls,tG_ls,tbnew,tC] = inversion_leastSquare_sparse_2019(u,v,z,dz,[0 0]);

    O_ls=[O_ls; tO_ls];
    G_ls=[G_ls; tG_ls];
    bnew=[bnew; tbnew'];
    C   =[C; tC];    
    
        % Display progress through the loop
    if jj==floor(length(seg_start)*.25)
        disp('25% of the way through separating ocean and glider velocity')
    elseif jj==floor(length(seg_start)*.5)
        disp('50% of the way through separating ocean and glider velocity')
    elseif jj==floor(length(seg_start)*.75)
        disp('75% of the way through separating ocean and glider velocity')
    end
    end %end segment loop
    clear u v z seg_time_start seg_time_end seg_count seg_time_ind tC tG_ls tO_ls tbnew

    save([data_path deployment '_ocean_glider_velo.mat'], 'O_ls', 'G_ls', 'bnew', 'C');
else

    load([data_path deployment '_ocean_glider_velo.mat'])

end

%% Split the ocean velocity profiles into E/W and N/S components
% u=real(O_ls) to get real (E/W) component
% v=imag(O_ls) to get imag (N/S) component
u=real(O_ls);
v=imag(O_ls);

%% Build grid using bnew and populate with ocean velocity
% Find min and max depth range for the grid from bnew
d_min=min(bnew);
d_max=max(bnew);

% Find the necessary number of depth bins
bin_num=((d_max-d_min)./dz)+1; % Add 1 to account for d_min bin
% Number of profile columns needed to fill with ocean velocity data
prof_num=length(seg_start);

% Create the grid
ugrid=nan([bin_num prof_num]);
vgrid=nan([bin_num prof_num]);
% Create a reference array for bin depth index
grid_bin=(d_min:dz:d_max)';

% Populate the grid
profile_ind=find(bnew==d_min); % Index of the start of each profile

for ii=1:length(profile_ind)-1
    % Grabs one profile of data
    prof_bins=bnew(profile_ind(ii):profile_ind(ii+1)-1)';
    prof_u   =u(profile_ind(ii):profile_ind(ii+1)-1)';
    prof_v   =v(profile_ind(ii):profile_ind(ii+1)-1)';
    
    % Finds the depths to at which to insert this data into the grid
    [~, grid_start]=min(prof_bins);
    [~, grid_end]=max(prof_bins);
    
    ugrid(grid_start:grid_end,ii)=prof_u;
    vgrid(grid_start:grid_end,ii)=prof_v;
end
    clear prof_u prof_v prof_bins grid_start grid_end

%% Plot the gridded ocean velocity data

% Get some colorbar limits. Define limits the same for u and v
tolerance=3; % standard deviation tolerance for the colorbar limits
velo=[u;v];
mvelo=nanmean(velo);
stdvelo=nanstd(velo);
clim=[mvelo-tolerance.*stdvelo mvelo+tolerance.*stdvelo];

% Set plotting constants
cmap=cmocean('delta');

% Plotting
figure(1)
pcolorjw(1:prof_num,-grid_bin,ugrid)
ax=gca;
hold on
cb=colorbar;
cb_label=get(cb,'Label');
set(cb_label,'String','East (+)/ West(-) Ocean velocity [m/s]');
ax.CLim=clim;
xlabel('Segment number')
ylabel('Depth [m]')
colormap(cmap)
title([deployment ': East (+)/ West(-) Ocean velocity [m/s]'],'Interpreter','none')
print([fig_path 'u_ocean_velocity_cs'],'-dpng')
close all
clear ax cb cb_label

figure(2)
pcolorjw(1:prof_num,-grid_bin,vgrid)
ax=gca;
hold on
cb=colorbar;
cb_label=get(cb,'Label');
set(cb_label,'String','North (+)/ South(-) Ocean velocity [m/s]');
ax.CLim=clim;
xlabel('Segment number')
ylabel('Depth [m]')
colormap(cmap)
title([deployment ': North (+)/ South(-) Ocean velocity [m/s]'],'Interpreter','none')
print([fig_path 'v_ocean_velocity_cs'],'-dpng')
close all
clear ax cb cb_label
