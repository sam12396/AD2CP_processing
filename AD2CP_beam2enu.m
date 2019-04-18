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
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' mission '/'];


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
rng_cells=ncread(data_path,'Average_Range',[1,1],[1,15])';

% Get the time variable
avg_time=ncread(data_path,'Average_MatlabTimeStamp');

% Get the beam velocity data
v1=double(ncread(data_path,'Average_VelBeam1')); % Forward X beam
v2=double(ncread(data_path,'Average_VelBeam2'));
v3=double(ncread(data_path,'Average_VelBeam3')); % Rear X beam
v4=double(ncread(data_path,'Average_VelBeam4'));

% Get magnetically corrected behavior data
head=double(ncread(data_path,'Average_Heading_Mag_Corr'));
pitch=double(ncread(data_path,'Average_Pitch'));
roll=double(ncread(data_path,'Average_Roll'));
pressure=double(ncread(data_path,'Average_Pressure'));
ensemble=double(ncread(data_path,'Average_EnsembleCount'));

%% Beam specific QA/QC
% Some quality control methods require beam velocities and thus, they need
% to be done before the conversion to ENU velocities. We do not need to
% worry about Correlation values because we are using the average data set

%% Beam Amplitude QA/QC
% Remove any bin and bin further away from it that returns an amplitude
% within 3dB of the noise floor. We also need to calculate the noise floor.
% This method is recommended by the manual
% Todd used a different method. Use SNR to do a similar check but say any
% bin with SNR<20 is bad
ampbeam1=double(ncread(data_path,'Average_AmpBeam1'));
ampbeam2=double(ncread(data_path,'Average_AmpBeam2'));
ampbeam3=double(ncread(data_path,'Average_AmpBeam3'));
ampbeam4=double(ncread(data_path,'Average_AmpBeam4'));

% Define the noise floor
nf=25; %dB. Reported by Nortek by sampling in air 

% SNR_dB=S_dB-N_dB use 13dB as theshold (Todd 2017). 
snr_threshold=13;
for ii=1:size(ampbeam1,1)
    if pitch(ii)>15 % Upcast so use beams 234
        % Calculate SNR for each beam
        snr2=(ampbeam2(ii,:)-nf);
        snr3=(ampbeam3(ii,:)-nf);
        snr4=(ampbeam4(ii,:)-nf);
        
        % If the value of SNR < threshold in a bin in any beam, make that
        % bin nan in all beams needed for this cast
        snr_ind=snr2<snr_threshold | snr3<snr_threshold | snr4<snr_threshold;
        ampbeam2(ii,snr_ind)=NaN;
        ampbeam3(ii,snr_ind)=NaN;
        ampbeam4(ii,snr_ind)=NaN;
    elseif pitch(ii)<-15 % Downcast so use beams 124
        % Calculate SNR for each beam
        snr1=(ampbeam1(ii,:)-nf);
        snr2=(ampbeam2(ii,:)-nf); 
        snr4=(ampbeam4(ii,:)-nf);
        
        % If the value of SNR < threshold in a bin in any beam, make that
        % bin nan in all beams needed for this cast        
        snr_ind=snr1<snr_threshold | snr2<snr_threshold | snr4<snr_threshold;
        ampbeam1(ii,snr_ind)=NaN;
        ampbeam2(ii,snr_ind)=NaN;
        ampbeam4(ii,snr_ind)=NaN;
    end
    if ii==floor(size(ampbeam1,1)*0.25)
        disp('25% of the way through SNR check')
    elseif ii==floor(size(ampbeam1,1)*0.50)
        disp('50% of the way through SNR check')
    elseif ii==floor(size(ampbeam1,1)*0.75)
        disp('75% of the way through SNR check')
    end
    clear snr_ind snr1 snr2 snr3 snr4
end

% Check how much data we lost in each beam
a=[];
for ii=1:size(ampbeam2,2)
    a(1,ii)=sum(isnan(ampbeam1(:,ii)))./length(ampbeam1).*100;
    a(2,ii)=sum(isnan(ampbeam2(:,ii)))./length(ampbeam1).*100;
    a(3,ii)=sum(isnan(ampbeam3(:,ii)))./length(ampbeam1).*100;
    a(4,ii)=sum(isnan(ampbeam4(:,ii)))./length(ampbeam1).*100;
end
figure(1)
hold on
scatter(1:15,a(1,:),'r');plot(1:15,a(1,:),'r');
scatter(1:15,a(2,:),'g');plot(1:15,a(2,:),'g');
scatter(1:15,a(3,:),'k');plot(1:15,a(3,:),'k');
scatter(1:15,a(4,:),'b');plot(1:15,a(4,:),'b');
title(['# of Nans by bins across all ensembles from SNR threshold of: ' num2str(snr_threshold) 'dB and Noise Floor of: ' num2str(nf) 'dB'])
legend('Beam1','Beam1','Beam2','Beam2','Beam3','Beam3','Beam4','Beam4','Location','NW')
ylabel('% NaN in data')
xlabel('Bin number')
print([fig_path 'nans_from_snr_threshold_' num2str(snr_threshold) '_dB'],'-dpng')
clear a 
%close all

% Remove extreme amplitude returns
% Set extreme threshold
inter_level=75;
for ii=1:size(ampbeam1,1)
    if pitch(ii)>15 % Upcast so use beams 234
        % If the return is above the threshold, flag the data
        amp2_ind= ampbeam2(ii,:)>inter_level;
        amp3_ind= ampbeam3(ii,:)>inter_level;
        amp4_ind= ampbeam4(ii,:)>inter_level;
        amp_ind=amp2_ind+amp3_ind+amp4_ind;
        ampbeam2(amp_ind~=0,:)=NaN;
        ampbeam3(amp_ind~=0,:)=NaN;
        ampbeam4(amp_ind~=0,:)=NaN;
        
    elseif pitch(ii)<-15 % Downcast so use beams 124
        % If the return is above the threshold, flag the data
        amp1_ind= ampbeam1(ii,:)>inter_level;
        amp2_ind= ampbeam2(ii,:)>inter_level;
        amp4_ind= ampbeam4(ii,:)>inter_level;
        amp_ind=amp1_ind+amp2_ind+amp4_ind;
        ampbeam1(amp_ind~=0,:)=NaN;
        ampbeam2(amp_ind~=0,:)=NaN;
        ampbeam4(amp_ind~=0,:)=NaN;
    end
    if ii==floor(size(ampbeam1,1)*0.25)
        disp('25% of the way through SNR check')
    elseif ii==floor(size(ampbeam1,1)*0.50)
        disp('50% of the way through SNR check')
    elseif ii==floor(size(ampbeam1,1)*0.75)
        disp('75% of the way through SNR check')
    end
    clear amp_ind amp1_ind amp2_ind amp3_ind amp4_ind
end
%% Apply qa/qc to actual beam velocities

qaqc=input('Perform QA/QC on the velocity data? Y=1 N=0');

if qaqc==1
    qaqc1_ind=isnan(ampbeam1);
    qaqc2_ind=isnan(ampbeam2);
    qaqc3_ind=isnan(ampbeam3);
    qaqc4_ind=isnan(ampbeam4);

    v1(qaqc1_ind)=NaN;
    v2(qaqc2_ind)=NaN;
    v3(qaqc3_ind)=NaN;
    v4(qaqc4_ind)=NaN;
    disp('QA/QC complete')
else
    disp('QA/QC skipped')
end

%% Convert beam velocity to ENU velocity

% Do pitch roll if statements outside of the for loop and trim data that
% way the loop doesnt take as long. Only take pitchs that indicate they
% are during a dive (>+-15). Only take points with roll inside +-10.
trim_ind= abs(pitch)>15 & abs(roll)<10;
trim_pitch=pitch(trim_ind);
trim_roll=roll(trim_ind);
trim_head=head(trim_ind);
trim_avg_time=avg_time(trim_ind);
trim_v1=v1(trim_ind,:);
trim_v2=v2(trim_ind,:);
trim_v3=v3(trim_ind,:);
trim_v4=v4(trim_ind,:);
trim_pressure=pressure(trim_ind);

kk=1;
for k=1:length(trim_avg_time) %loop through each ensemble
    
V1=trim_v1(k,:);
V2=trim_v2(k,:);
V3=trim_v3(k,:);
V4=trim_v4(k,:);

% Concat all up/down cast beam data
beam124 = [V1;V2;V4];
beam234 = [V2;V3;V4];
    
% Calculate the vertical displacement of the bins adjusting for pitch and 
% roll
zrng = cell_vert(trim_pitch(k),trim_roll(k),rng_cells); 
zrng = trim_pressure(k) + zrng; % this is a real bin depth now

% Convert the beam velocity data to ENU velocity data.
% Only convert the front beams if downcast and back beams for upcast
    if trim_pitch(k) < -15
        [enu124] = Xform(beam124,T124,trim_head(k),trim_pitch(k),trim_roll(k),1);
        X.r(kk,:) = zrng.';  % Creates a field with the real bin depth
        X.e(kk,:)=enu124(1,:);
        X.n(kk,:)=enu124(2,:);
        X.u(kk,:)=enu124(3,:);        
    elseif trim_pitch(k) > 15
        [enu234] = Xform(beam234,T234,trim_head(k),trim_pitch(k),trim_roll(k),1);
        X.r(kk,:) = zrng.';  % Creates a field with the real bin depth
        X.e(kk,:)=enu234(1,:);
        X.n(kk,:)=enu234(2,:);
        X.u(kk,:)=enu234(3,:);        
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

% Store trimmed time variable
X.t=trim_avg_time;

%% Write the converted data to the file
ens_num=length(trim_avg_time);
bin_num=size(trim_v1,2);

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
        {'Trimmed_Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_NVelocity','Dimensions',...
        {'Trimmed_Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_UVelocity','Dimensions',...
        {'Trimmed_Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_Depth','Dimensions',...
        {'Trimmed_Ensemble_number',ens_num,'Bins',bin_num},'Datatype','double');
    nccreate(data_path,'Average_TrimMatlabTimeStamp','Dimensions',...
        {'Trimmed_Ensemble_number',ens_num},'Datatype','double');    
end

ncwrite(data_path,'Average_EVelocity',X.e);
ncwrite(data_path,'Average_NVelocity',X.n);
ncwrite(data_path,'Average_UVelocity',X.u);
ncwrite(data_path,'Average_Depth',X.r);
ncwrite(data_path,'Average_TrimMatlabTimeStamp',X.t);

% Add comments and units
ncwriteatt(data_path,'Average_EVelocity','Units','m/s');
ncwriteatt(data_path,'Average_NVelocity','Units','m/s');
ncwriteatt(data_path,'Average_UVelocity','Units','m/s');
ncwriteatt(data_path,'Average_Depth','Units','m')
ncwriteatt(data_path,'Average_TrimMatlabTimeStamp','Units','days');

ncwriteatt(data_path,'Average_EVelocity','Comments',['East/West velocity. Converted from beam velocity using Xform.m. SNR threshold of ' num2str(snr_threshold) ' was used']);
ncwriteatt(data_path,'Average_NVelocity','Comments',['North/South velocity. Converted from beam velocity using Xform.m. SNR threshold of ' num2str(snr_threshold) ' was used']);
ncwriteatt(data_path,'Average_UVelocity','Comments',['Up/Down velocity. Converted from beam velocity using Xform.m. SNR threshold of ' num2str(snr_threshold) ' was used']);
ncwriteatt(data_path,'Average_Depth','Comments','Actual depth of bottom of each bin. From AD2CP_beam2enu.m.')
ncwriteatt(data_path,'Average_TrimMatlabTimeStamp','Comments','Matlab timestamps of profiles remaining after pitch/roll trimming');

