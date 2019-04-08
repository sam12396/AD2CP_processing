%% Calculate Richardson number from AD2CP and glider data
% By using the shear data calculcated from data collected by the AD2CP and
% buoyancy frequency calculated from glider data, richardson number can be
% found. The AD2CP data is split by segments of the glider mission and
% depth bins. So, the bouyancy frequency variable will need averaged in the
% same way so that shear and buoyancy fequency have the same shape.

% Sam Coakley
% 2/28/2019

%% Starters
clear all;close all;

deployment='ru31_467_2016_01_13';
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' deployment '/'];
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment '/'];

a_data=load([data_path deployment '_ocean_shear.mat']);
g_data=load([data_path deployment '_dgroup.mat']);

%% Pull out AD2CP vars

% Depth of tops of bins
grid_bin=a_data.grid_bin;
dz=grid_bin(2)-grid_bin(1);
% Shear magnitude grid
s_maggrid=a_data.s_maggrid;

%% Calculate buoyancy frequency N2

% Calculate profile by profile N2 not segment by segment
% Sorts the variables so depth is monotonic increasing in each profile
% Removes Nans from glider data
[g_time, g_dep, N2]=calc_dgroup_N2(g_data.dgroup);

%% Grid N2 data [bin x profile]
% Need to make N2 a segment variable by averaging together each profile in
% a segment at the depth of the bin centers
% Start by putting the glider data on a bin by profile grid

% Find total number of profiles
prof_count=sum(diff(g_dep)<0)+1;
% Write number of bins
bin_num=length(grid_bin);
% Create grid to fill with N2 bin averaged data
N2_grid=nan(bin_num,prof_count);
profile_time=nan(1,prof_count);
% Make profile index of when profiles start
profile_ind=[0; find(diff(g_dep)<0)]+1; 

% Populate N2 grid 
for ii=1:prof_count-1
    % Grab one profile of data
    prof_dep=g_dep(profile_ind(ii):(profile_ind(ii+1)-1));
    prof_N2 =N2(profile_ind(ii):(profile_ind(ii+1)-1));
    prof_time=g_time(profile_ind(ii):(profile_ind(ii+1)-1));
    
    % Average data together that is inside the same depth bin
    for jj=1:bin_num
        % Find the indexs of data inside one bin. grid_bin is bin centers
        % so we find values around it and average into the bin
        bin_ind=find(prof_dep>=(grid_bin(jj)-(dz./2)) & prof_dep<(grid_bin(jj))+(dz./2));
        % Average data in that bin
        N2_grid(jj,ii)=nanmean(prof_N2(bin_ind));  
        % N2_grid may not be filled to the last bin because it is
        % calculated from the ctd on the glider and the AD2CP samples below
        % the max depth of the glider
    end    
    profile_time(1,ii)=nanmean(prof_time);
    clear prof_*
end

%% Grid N2 now by [bin x segment average]

% Get start and end times for each glider segment
segs=g_data.dgroup.startDatenums;
sege=g_data.dgroup.endDatenums;

% Create grid to fill with N2 segement averaged data
N2_seggrid=nan(bin_num,length(segs));

% Average together all profiles in a segment
for jj=1:length(segs)
    % Grab one segment of data
    seg_time=profile_time>= segs(jj) & profile_time< sege(jj);
    % Average the segment and populate the grid
    N2_seggrid(:,jj) =nanmean(N2_grid(:,seg_time),2);
end

%% Build grid of profile N2 with full depth resolution
% Create 1m resolution bin centers
bin1m=0.5:1:ceil(max(g_dep))-0.5;

% Preallocate memory
N2_grid1m=nan(length(bin1m),length(profile_ind));
for ii=1:length(profile_ind)-1
    % Grab one profile of data
    prof_dep=g_dep(profile_ind(ii):(profile_ind(ii+1)-1));
    prof_N2 =N2(profile_ind(ii):(profile_ind(ii+1)-1));
    prof_time=g_time(profile_ind(ii):(profile_ind(ii+1)-1));
    
    % Average data together that is inside the same depth bin
    for jj=1:length(bin1m)
        % Find the indexs of data inside one bin. bin1m is bin centers
        % so we find values around it and average into the bin
        bin_ind=find(prof_dep>=(bin1m(jj)-(1./2)) & prof_dep<(bin1m(jj))+(1./2));
        % Average data in that bin
        N2_grid1m(jj,ii)=nanmean(prof_N2(bin_ind));  
        % N2_grid1m may not be filled to the last bin because it is
        % calculated from the ctd on the glider and the AD2CP samples below
        % the max depth of the glider. The calculation of N2 also chops one
        % depth off each profile
    end    
    profile_time(ii)=nanmean(prof_time);
    clear prof_*
end

%% Calculate Richardson number
% Ri=N^2/(dU/dz)^2
% Reduced Shear= 4N^2 - (dU/dz)^2
ri_grid=N2_seggrid./(s_maggrid).^2;
segment_start=segs;
red_shear=4.*N2_seggrid - s_maggrid.^2;
save([data_path deployment '_rich_num_grid.mat'], 'red_shear', 'ri_grid', 'N2_grid', 'profile_time','segment_start', 'N2_grid1m', 'bin1m')