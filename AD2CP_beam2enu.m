%% Conversion of beam velocity to ENU velocities for AD2CP

%% Notes
% This script is for netcdf AD2CP data that has already had the magnetic
% corrections made to the heading data. This processing script should only
% be used after the data was processed by the AD2CP_mag_correct.m script.
% This script will convert the beam velocities in the data to East North Up
% velocities (ENU). The ENU velocities will then be added as a new variable
% to the netcdf.

% This is step 2 in AD2CP data processing. The velocity data still
% represents the sum of the ocean velocity observed and the glider
% velocity.

% Sam Coakley
% 2/27/19

%% Setting constants and loading data
clear all; close all;
mission='ru31_467_2016_01_13';
fname=[mission '_AD2CP.nc'];
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' mission '/' fname];

% These are the transformation matricies for up and down cast. Not sure
% which is which so I need to find some documentation.
% Left, Rear, Right beams
T234=[ 0.50552979193390  -1.35634170490926   0.50552979193390;...
  -1.18310079157625                  0   1.18310079157625;...
   0.55168895948125                  0   0.55168895948125];

% Forward, left, right beams
T124 = [1.35634170490926  -0.50552979193390  -0.50552979193390;...
                  0  -1.18310079157625   1.18310079157625;...
                  0   0.55168895948125   0.55168895948125];

% Get the bottom of each depth bin
rng_cells=ncread(data_path,'Average_Range',[1,1],[1,15],[1,1])';

% Get the time variable
avg_time=ncread(data_path,'Average_MatlabTimeStamp');

% Get the beam velocity data
v1=double(ncread(data_path,'Average_VelBeam1')); % Forward X beam
v2=double(ncread(data_path,'Average_VelBeam2'));
v3=double(ncread(data_path,'Average_VelBeam3')); % Rear X beam
v4=double(ncread(data_path,'Average_VelBeam4'));

% Get magnetically corrected behavior data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do I need a mag corr for pitch, roll ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
head=double(ncread(data_path,'Average_Heading_Mag_Corr'));
pitch=double(ncread(data_path,'Average_Pitch'));
roll=double(ncread(data_path,'Average_Roll'));
pressure=double(ncread(data_path,'Average_Pressure'));
ensemble=double(ncread(data_path,'Average_EnsembleCount'));

%% Convert beam velocity to ENU velocity

kk=1;
for k=1:length(avg_time) %loop through each ensemble
    
V1=v1(k,:);
V2=v2(k,:);
V3=v3(k,:);
V4=v4(k,:);

% Concat all up/down cast beam data
beam124 = [V1;V2;V4];
beam234 = [V2;V3;V4];

% Convert the beam velocity data to ENU velocity data
% Is the setting of 1 for the last argument correct? I do not understand
% that parameter. Why was the bottom track 0
[enu124] = Xform(beam124,T124,head(k),pitch(k),roll(k),1);
[enu234] = Xform(beam234,T234,head(k),pitch(k),roll(k),1);

% Calculate the vertical displacement of the bins adjusting for pitch and 
% roll
zrng = cell_vert(pitch(k),roll(k),rng_cells); 
zrng = pressure(k) + zrng; % this is a real bin depth now

            % Creates a field with the real bin depth
            X.r(kk,:) = zrng.'; 
            % Change the beam data being used for up or down cast. The beam
            % orientation is 1=Forward; 2 and 4=Port and Starboard; 3= Aft.
            % +-15° removes data from when the glider is inflecting. Skip
            % anything with too great of a roll.
            if pitch(k) < -15 && abs(roll(k))>10
            X.e(kk,:)=enu124(1,:);
            X.n(kk,:)=enu124(2,:);
            X.u(kk,:)=enu124(3,:);
            elseif pitch(k) > 15 && abs(roll(k))>10
            X.e(kk,:)=enu234(1,:);
            X.n(kk,:)=enu234(2,:);
            X.u(kk,:)=enu234(3,:);
            else
            X.e(kk,:)=nan(size(enu124));
            X.n(kk,:)=nan(size(enu124));
            X.u(kk,:)=nan(size(enu124));                
            end
kk=kk+1;

    % Display progress through the loop
    if k==floor(length(avg_time)*.25)
        disp('25% of the way through beam2enu conversion')
    elseif k==floor(length(avg_time)*.5)
        disp('50% of the way through beam2enu conversion')
    elseif k==floor(length(avg_time)*.75)
        disp('75% of the way through beam2enu conversion')
    end
end %end ensemble loop

%% Write the converted data to the file
ens_num=length(avg_time);
bin_num=size(v1,2);

% Stop variables from being created if this script has been run before and
% the variables already exist in the file.
info=ncinfo(data_path);
for ii=1:length(info.Variables)
    % Check to see if the variable appears in the data file.
    % Returns a 1 if the strings are equal ie the variable already exists.
    create_vars=strcmp(info.Variables(ii).Name,'Average_EVelocity');
    
    % Leave the loop if the variable already exists
    if create_vars==1
        break
    end
end

% Skip creation of the variables if they already exist in the file
if create_vars~=1
    % Create/write the data from the structure to the .nc file
    nccreate(data_path,'Average_EVelocity','Dimensions',...
        {'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_NVelocity','Dimensions',...
        {'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_UVelocity','Dimensions',...
        {'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_Depth','Dimensions',...
        {'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
end

ncwrite(data_path,'Average_EVelocity',X.e);
ncwrite(data_path,'Average_NVelocity',X.n);
ncwrite(data_path,'Average_UVelocity',X.u);
ncwrite(data_path,'Average_Depth',X.r);

% Add comments and units
ncwriteatt(data_path,'Average_EVelocity','Units','m/s');
ncwriteatt(data_path,'Average_NVelocity','Units','m/s');
ncwriteatt(data_path,'Average_UVelocity','Units','m/s');
ncwriteatt(data_path,'Average_Depth','Units','m')

ncwriteatt(data_path,'Average_EVelocity','Comments','East/West velocity. Converted from beam velocity using Xform.m');
ncwriteatt(data_path,'Average_NVelocity','Comments','North/South velocity. Converted from beam velocity using Xform.m');
ncwriteatt(data_path,'Average_UVelocity','Comments','Up/Down velocity. Converted from beam velocity using Xform.m');
ncwriteatt(data_path,'Average_Depth','Comments','Actual depth of bottom of each bin. From AD2CP_beam2enu.m.')

