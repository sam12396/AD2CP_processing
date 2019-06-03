%% Magnetic Correction processing for AD2CP
% This script is used to concatenate all of the data files into one and add
% a new field of magnetically-corrected heading to the data set.
% This script is step 1 in processing glider mounted AD2CP data.

% 2/13/19
% Sam Coakley

%%
clear all;close all;

% Path to directory with ADCP data.
deployment='ru31_467_2016_01_13/';
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment];
g_fname='ru31_467_dgroup_160113.mat';
% glider_data=[data_path 'ru31_467_dgroup_160113.mat'];
%% Concatenate all files into one structure

% Load all the separate data files. They are likely just separated for size
% reasons because there is no indication otherwise.
a1=load([data_path 'Glider.ad2cp.00000_1.mat']);
a2=load([data_path 'Glider.ad2cp.00000_2.mat']);
a3=load([data_path 'Glider.ad2cp.00000_3.mat']);
a4=load([data_path 'Glider.ad2cp.00000_4.mat']);
a5=load([data_path 'Glider.ad2cp.00000_5.mat']);
a6=load([data_path 'Glider.ad2cp.00000_6.mat']);
a7=load([data_path 'Glider.ad2cp.00000_7.mat']);
a8=load([data_path 'Glider.ad2cp.00000_8.mat']);
a9=load([data_path 'Glider.ad2cp.00000_9.mat']);

% Concatenate the data into one large matrix
at=[a1 a2 a3 a4 a5 a6 a7 a8 a9];
clear a1 a2 a3 a4 a5 a6 a7 a8 a9

%% Concatenate data from all files into single variables

headt=[];
pitcht=[]; 
rollt=[]; 
pressuret=[]; 
xt=[]; 
yt=[]; 
zt=[]; 
timet=[];
for k=1:9 % Loop through each data file in the structure
    
    % Time and behavior data
    time=double(at(k).Data.Average_MatlabTimeStamp);
    head=double(at(k).Data.Average_Heading); % AD2CP compass heading (may need hard iron correction)
    pitch=double(at(k).Data.Average_Pitch); % AD2CP pitch (may need hard iron correction)
    roll=double(at(k).Data.Average_Roll); % AD2CP roll (may need hard iron correction)
    pressure=double(at(k).Data.Average_Pressure); % AD2CP pressure sensor

    % Magnetometer data
    x=double(at(k).Data.Average_MagnetometerX);
    y=double(at(k).Data.Average_MagnetometerY);
    z=double(at(k).Data.Average_MagnetometerZ);

    % Store in temporary variables
    headt=[headt head']; 
    pitcht=[pitcht pitch']; 
    rollt=[rollt roll']; 
    pressuret=[pressuret pressure']; 
    xt=[xt x']; 
    yt=[yt y']; 
    zt=[zt z']; 
    timet=[timet time'];
end % k loop through all files

% Rename and transpose the temporary variables
head=headt'; 
pitch=pitcht'; 
roll=rollt'; 
pressure=pressuret'; 
x=xt'; 
y=yt'; 
z=zt';
time=timet';

clear headt pitcht rollt pressuret xt yt zt timet timeBT

XYZ_original=[x y z];
%% Pitch Dependent magnetic heading correction
% The pitch of the glider changes by moving the battery pack. The battery
% pack produces a magnetic field so when it moves, the field moves and this
% effects our compass data. Having all of the data from a mission loaded in
% for this correction ensures the best fit possible.

pitch_ranges=(-50:10:50);

for k=1:length(pitch_ranges)-1

    % Find where pitch falls into this particular range    
    ii=pitch >= pitch_ranges(k) & pitch < pitch_ranges(k+1);

    % Index the magnetic field variables to only be referencing the current
    % pitch range
    XYZ1=[x(ii) y(ii) z(ii)];
 
    % Fits the given cartesian coordinates to an ellipsoid and gives the
    % center of that shape as the first output.
    [offset, radii, evecs, evals, pars ] = ellipsoid_fit_tnm( XYZ1 );
 
    % Creates xyz variables that are offset by the center of the
    % magnetic field recorded in the data
    x1 = x(ii)-offset(1);
    y1 = y(ii)-offset(2);
    z1 = z(ii)-offset(3);

    % Find the center of the new magnetic data that should be closer 0,0,z
    % than the original data.
    [new_center, radii, evecs, evals, pars ] = ellipsoid_fit_tnm([x1 y1 z1]);
    
    % If the center of the offset magnetic data is not near 0,0 in the x
    % and y planes, we keep the original values.
    if abs(new_center(1)) > 150 || abs(new_center(2) > 150)
        x1=x(ii); y1=y(ii);z1=z(ii);
    end

    % Replace the indexed values with the new corrected values
    % If the previous if statement is true, nothing happens here, x(ii)=x1
    % which = x(ii).
    % If it is not, then x(ii)=x(ii)-offset1 etc.
    x(ii)=x1; 
    y(ii)=y1; 
    z(ii)=z1;
    clear x1 y1 z1 offset radii evecs evals pars
end

XYZ_final=[x y z];

%% Re-calculate heading with magnetometer corrections
for k=1:length(pressure)
% Final heading
headingf(k) = CalcMidlifeHeading(XYZ_final(k,:), pitch(k), roll(k), 1);
% Original heading
headingo(k) = CalcMidlifeHeading(XYZ_original(k,:), pitch(k), roll(k), 1);
end

% Plots the corrected (blue) and uncorrected (red) heading values to check
% the correction was successful.
close all;plot(XYZ_original(:,1),XYZ_original(:,2),'.r',x,y,'.b')
corrected_heading = headingf';
%% Put correction back into individual data structures

% Create new data structure
at2=at;

% Find the total number of Average_Heading data points
for k=1:9
    lf=length(at2(k).Data.Average_Heading);
    tfl(k)=lf;
end
    % Create a indexing variable for the corrected heading
    cumt=cumsum(tfl);

    % The full length of all the data arrays
    total_len=sum(tfl);
% Add a field with the corrected heading data
for k=1:9
   if k==1
    at2(k).Data.Average_Heading_Mag_Corr = single(corrected_heading(1:cumt(k)));
   else
    at2(k).Data.Average_Heading_Mag_Corr = single(corrected_heading(cumt(k-1)+1:cumt(k)));
   end
end

clear at
%% Concatenate all variables from different files into full data sets
% Going to make arrays of each field with all the mission data. The end
% goal is to put these into a netcdf so that we can call segments easily
% when loading.

bin_num=size(at2(1).Data.Average_VelBeam1,2);
% Preallocate memory for each variable
Average_VelBeam1=nan(total_len, bin_num,'single');
Average_VelBeam2=nan(total_len, bin_num,'single');
Average_VelBeam3=nan(total_len, bin_num,'single');
Average_VelBeam4=nan(total_len, bin_num,'single');
Average_CorBeam1=nan(total_len, bin_num,'single');
Average_CorBeam2=nan(total_len, bin_num,'single');
Average_CorBeam3=nan(total_len, bin_num,'single');
Average_CorBeam4=nan(total_len, bin_num,'single');
Average_AmpBeam1=nan(total_len, bin_num,'single');
Average_AmpBeam2=nan(total_len, bin_num,'single');
Average_AmpBeam3=nan(total_len, bin_num,'single');
Average_AmpBeam4=nan(total_len, bin_num,'single');
Average_TimeStamp=nan(total_len,1,'double');
Average_MatlabTimeStamp=nan(total_len,1,'double');
Average_SerialNumber=nan(total_len,1,'single');
Average_SpeedOfSound=nan(total_len,1,'single');
Average_Error=nan(total_len,1,'single');
Average_Status=nan(total_len,1,'single');
Average_CellSize=nan(total_len,1,'single');
Average_Blanking=nan(total_len,1,'single');
Average_NominalCor=nan(total_len,1,'single');
Average_Battery=nan(total_len,1,'single');
Average_MagnetometerX=nan(total_len,1,'single');
Average_MagnetometerY=nan(total_len,1,'single');
Average_MagnetometerZ=nan(total_len,1,'single');
Average_AccelerometerX=nan(total_len,1,'single');
Average_AccelerometerY=nan(total_len,1,'single');
Average_AccelerometerZ=nan(total_len,1,'single');
Average_AmbiguityVel=nan(total_len,1,'single');
Average_TransmitEnergy=nan(total_len,1,'single');
Average_PowerLevel=nan(total_len,1,'single');
Average_Physicalbeam=nan(total_len,4,'single');
Average_MagnetometerTemperature=nan(total_len,1,'single');
Average_RTCTemperature=nan(total_len,1,'single');
Average_PressureSensorTemperature=nan(total_len,1,'single');
Average_EnsembleCount=nan(total_len,1,'single');
Average_WaterTemperature=nan(total_len,1,'single');
Average_Pressure=nan(total_len,1,'single');
Average_Heading=nan(total_len,1,'single');
Average_Pitch=nan(total_len,1,'single');
Average_Roll=nan(total_len,1,'single');
Average_Range=nan(length(at2),bin_num,'single');
Average_Heading_Mag_Corr=nan(total_len,1,'single');

Units=at2(1).Data.Units;
Comments=at2(1).Data.Comments;


% Loop through each file 

for ii=1:length(at2)
    % Fill variables with data from structure
    if ii==1
        Average_VelBeam1(1:cumt(ii),:)=at2(ii).Data.Average_VelBeam1;
        Average_VelBeam2(1:cumt(ii),:)=at2(ii).Data.Average_VelBeam2;
        Average_VelBeam3(1:cumt(ii),:)=at2(ii).Data.Average_VelBeam3;
        Average_VelBeam4(1:cumt(ii),:)=at2(ii).Data.Average_VelBeam4;
        Average_CorBeam1(1:cumt(ii),:)=at2(ii).Data.Average_CorBeam1;
        Average_CorBeam2(1:cumt(ii),:)=at2(ii).Data.Average_CorBeam2;
        Average_CorBeam3(1:cumt(ii),:)=at2(ii).Data.Average_CorBeam3;
        Average_CorBeam4(1:cumt(ii),:)=at2(ii).Data.Average_CorBeam4;
        Average_AmpBeam1(1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam1;
        Average_AmpBeam2(1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam2;
        Average_AmpBeam3(1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam3;
        Average_AmpBeam4(1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam4;
        Average_TimeStamp(1:cumt(ii),:)=at2(ii).Data.Average_TimeStamp;
        Average_MatlabTimeStamp(1:cumt(ii),:)=at2(ii).Data.Average_MatlabTimeStamp;
        Average_SerialNumber(1:cumt(ii),:)=at2(ii).Data.Average_SerialNumber;
        Average_SpeedOfSound(1:cumt(ii),:)=at2(ii).Data.Average_SpeedOfSound;
        Average_Error(1:cumt(ii),:)=at2(ii).Data.Average_Error;
        Average_Status(1:cumt(ii),:)=at2(ii).Data.Average_Status;
        Average_CellSize(1:cumt(ii),:)=at2(ii).Data.Average_CellSize;
        Average_Blanking(1:cumt(ii),:)=at2(ii).Data.Average_Blanking;
        Average_NominalCor(1:cumt(ii),:)=at2(ii).Data.Average_NominalCor;
        Average_Battery(1:cumt(ii),:)=at2(ii).Data.Average_Battery;
        Average_MagnetometerX(1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerX;
        Average_MagnetometerY(1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerY;
        Average_MagnetometerZ(1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerZ;
        Average_AccelerometerX(1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerX;
        Average_AccelerometerY(1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerY;
        Average_AccelerometerZ(1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerZ;
        Average_AmbiguityVel(1:cumt(ii),:)=at2(ii).Data.Average_AmbiguityVel;
        Average_TransmitEnergy(1:cumt(ii),:)=at2(ii).Data.Average_TransmitEnergy;
        Average_PowerLevel(1:cumt(ii),:)=at2(ii).Data.Average_PowerLevel;
        Average_Physicalbeam(1:cumt(ii),:)=at2(ii).Data.Average_Physicalbeam;
        Average_MagnetometerTemperature(1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerTemperature;
        Average_RTCTemperature(1:cumt(ii),:)=at2(ii).Data.Average_RTCTemperature;
        Average_PressureSensorTemperature(1:cumt(ii),:)=at2(ii).Data.Average_PressureSensorTemperature;
        Average_EnsembleCount(1:cumt(ii),:)=at2(ii).Data.Average_EnsembleCount;
        Average_WaterTemperature(1:cumt(ii),:)=at2(ii).Data.Average_WaterTemperature;
        Average_Pressure(1:cumt(ii),:)=at2(ii).Data.Average_Pressure;
        Average_Heading(1:cumt(ii),:)=at2(ii).Data.Average_Heading;
        Average_Pitch(1:cumt(ii),:)=at2(ii).Data.Average_Pitch;
        Average_Roll(1:cumt(ii),:)=at2(ii).Data.Average_Roll;
        Average_Range(ii,:)=at2(ii).Data.Average_Range;
        Average_Heading_Mag_Corr(1:cumt(ii),:)=at2(ii).Data.Average_Heading_Mag_Corr;
    else
        Average_VelBeam1(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_VelBeam1;
        Average_VelBeam2(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_VelBeam2;
        Average_VelBeam3(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_VelBeam3;
        Average_VelBeam4(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_VelBeam4;
        Average_CorBeam1(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_CorBeam1;
        Average_CorBeam2(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_CorBeam2;
        Average_CorBeam3(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_CorBeam3;
        Average_CorBeam4(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_CorBeam4;
        Average_AmpBeam1(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam1;
        Average_AmpBeam2(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam2;
        Average_AmpBeam3(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam3;
        Average_AmpBeam4(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AmpBeam4;
        Average_TimeStamp(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_TimeStamp;
        Average_MatlabTimeStamp(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_MatlabTimeStamp;
        Average_SerialNumber(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_SerialNumber;
        Average_SpeedOfSound(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_SpeedOfSound;
        Average_Error(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Error;
        Average_Status(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Status;
        Average_CellSize(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_CellSize;
        Average_Blanking(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Blanking;
        Average_NominalCor(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_NominalCor;
        Average_Battery(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Battery;
        Average_MagnetometerX(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerX;
        Average_MagnetometerY(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerY;
        Average_MagnetometerZ(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerZ;
        Average_AccelerometerX(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerX;
        Average_AccelerometerY(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerY;
        Average_AccelerometerZ(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AccelerometerZ;
        Average_AmbiguityVel(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_AmbiguityVel;
        Average_TransmitEnergy(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_TransmitEnergy;
        Average_PowerLevel(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_PowerLevel;
        Average_Physicalbeam(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Physicalbeam;
        Average_MagnetometerTemperature(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_MagnetometerTemperature;
        Average_RTCTemperature(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_RTCTemperature;
        Average_PressureSensorTemperature(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_PressureSensorTemperature;
        Average_EnsembleCount(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_EnsembleCount;
        Average_WaterTemperature(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_WaterTemperature;
        Average_Pressure(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Pressure;
        Average_Heading(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Heading;
        Average_Pitch(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Pitch;
        Average_Roll(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Roll;
        Average_Range(ii,:)=at2(ii).Data.Average_Range;
        Average_Heading_Mag_Corr(cumt(ii-1)+1:cumt(ii),:)=at2(ii).Data.Average_Heading_Mag_Corr;
    end
end

%% Make a netcdf with all these variables

% Need to use one nccreate per variable
% nccreate('filename','variablename','Dimensions',{'dimname1',dimvalue1,'dimname2',dimvalue2})
% Define dimensions needed
ens_num=length(Average_AccelerometerX);
if ~isfile([data_path deployment(1:end-1) '_AD2CP.nc'])
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AccelerometerX','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_AccelerometerY','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AccelerometerZ','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AmbiguityVel','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AmpBeam1','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AmpBeam2','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AmpBeam3','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_AmpBeam4','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_Battery','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_Blanking','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_CellSize','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_CorBeam1','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_CorBeam2','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_CorBeam3','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_CorBeam4','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_EnsembleCount','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Error','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Heading','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_Heading_Mag_Corr','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_MagnetometerTemperature','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_MagnetometerX','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_MagnetometerY','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_MagnetometerZ','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_MatlabTimeStamp','Dimensions',{'Ensemble_number',ens_num},'Datatype','double');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_NominalCor','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Physicalbeam','Dimensions',{'Ensemble_number',ens_num,'Beam_num',4},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Pitch','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_PowerLevel','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Pressure','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_PressureSensorTemperature','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_RTCTemperature','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Range','Dimensions',{'Files',9,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_Roll','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_SerialNumber','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_SpeedOfSound','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...
    'Average_Status','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_TimeStamp','Dimensions',{'Ensemble_number',ens_num},'Datatype','double');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_TransmitEnergy','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_VelBeam1','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_VelBeam2','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_VelBeam3','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_VelBeam4','Dimensions',{'Ensemble_number',ens_num,'Bins',bin_num},'Datatype','single');
nccreate([data_path deployment(1:end-1) '_AD2CP.nc'], ...    
    'Average_WaterTemperature','Dimensions',{'Ensemble_number',ens_num},'Datatype','single');
end

%% Use this to write the data to the file
%ncwrite('filename','variablename',data)
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerX',Average_AccelerometerX);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerY',Average_AccelerometerY);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerZ',Average_AccelerometerZ);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmbiguityVel',Average_AmbiguityVel);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam1',Average_AmpBeam1);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam2',Average_AmpBeam2);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam3',Average_AmpBeam3);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam4',Average_AmpBeam4);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Battery',Average_Battery);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Blanking',Average_Blanking);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CellSize',Average_CellSize);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam1',Average_CorBeam1);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam2',Average_CorBeam2);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam3',Average_CorBeam3);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam4',Average_CorBeam4);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_EnsembleCount',Average_EnsembleCount);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Error',Average_Error);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Heading',Average_Heading);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Heading_Mag_Corr',Average_Heading_Mag_Corr);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerTemperature',Average_MagnetometerTemperature);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerX',Average_MagnetometerX);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerY',Average_MagnetometerY);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerZ',Average_MagnetometerZ);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MatlabTimeStamp',Average_MatlabTimeStamp);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_NominalCor',Average_NominalCor);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Physicalbeam',Average_Physicalbeam);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Pitch',Average_Pitch);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PowerLevel',Average_PowerLevel);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Pressure',Average_Pressure);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PressureSensorTemperature',Average_PressureSensorTemperature);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_RTCTemperature',Average_RTCTemperature);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Range',Average_Range);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Roll',Average_Roll);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_SerialNumber',Average_SerialNumber);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_SpeedOfSound',Average_SpeedOfSound);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Status',Average_Status);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_TimeStamp',Average_TimeStamp);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_TransmitEnergy',Average_TransmitEnergy);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam1',Average_VelBeam1);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam2',Average_VelBeam2);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam3',Average_VelBeam3);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam4',Average_VelBeam4);
ncwrite([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_WaterTemperature',Average_WaterTemperature);

%% Use this to write attributes to the file ie Units and Comments
%ncwriteatt(filename,location,attname,attvalue)

% Write the comments from the original data into the netcdf
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam1','Comments',Comments.Average_CorBeam1);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam2','Comments',Comments.Average_CorBeam2);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam3','Comments',Comments.Average_CorBeam3);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam4','Comments',Comments.Average_CorBeam4);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_TimeStamp','Comments',Comments.Average_TimeStamp);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MatlabTimeStamp','Comments',Comments.Average_MatlabTimeStamp);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_SerialNumber','Comments',Comments.Average_SerialNumber);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Error','Comments',Comments.Average_Error);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Status','Comments',Comments.Average_Status);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CellSize','Comments',Comments.Average_CellSize);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Blanking','Comments',Comments.Average_Blanking);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerX','Comments',Comments.Average_MagnetometerX);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerY','Comments',Comments.Average_MagnetometerY);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerZ','Comments',Comments.Average_MagnetometerZ);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerX','Comments',Comments.Average_AccelerometerX);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerY','Comments',Comments.Average_AccelerometerY);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AccelerometerZ','Comments',Comments.Average_AccelerometerZ);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmbiguityVel','Comments',Comments.Average_AmbiguityVel);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_TransmitEnergy','Comments',Comments.Average_TransmitEnergy);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PowerLevel','Comments',Comments.Average_PowerLevel);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Physicalbeam','Comments',Comments.Average_Physicalbeam);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_RTCTemperature','Comments',Comments.Average_RTCTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PressureSensorTemperature','Comments',Comments.Average_PressureSensorTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Heading_Mag_Corr','Comments','Corrected using methods in AD2CP_detail_process.m');


% Write the units from the original data into the netcdf
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam1','Units',Units.Average_VelBeam1);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam2','Units',Units.Average_VelBeam2);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam3','Units',Units.Average_VelBeam3);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_VelBeam4','Units',Units.Average_VelBeam4);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam1','Units',Units.Average_CorBeam1);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam2','Units',Units.Average_CorBeam2);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam3','Units',Units.Average_CorBeam3);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CorBeam4','Units',Units.Average_CorBeam4);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam1','Units',Units.Average_AmpBeam1);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam2','Units',Units.Average_AmpBeam2);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam3','Units',Units.Average_AmpBeam3);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmpBeam4','Units',Units.Average_AmpBeam4);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_TimeStamp','Units',Units.Average_TimeStamp);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MatlabTimeStamp','Units',Units.Average_MatlabTimeStamp);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_SpeedOfSound','Units',Units.Average_SpeedOfSound);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_CellSize','Units',Units.Average_CellSize);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Blanking','Units',Units.Average_Blanking);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_NominalCor','Units',Units.Average_NominalCor);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Battery','Units',Units.Average_Battery);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_AmbiguityVel','Units',Units.Average_AmbiguityVel);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PowerLevel','Units',Units.Average_PowerLevel);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Pitch','Units',Units.Average_Pitch);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Pressure','Units',Units.Average_Pressure);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Heading','Units',Units.Average_Heading);
% Same units as heading
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Heading_Mag_Corr','Units',Units.Average_Heading);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_MagnetometerTemperature','Units',Units.Average_MagnetometerTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_RTCTemperature','Units',Units.Average_RTCTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_PressureSensorTemperature','Units',Units.Average_PressureSensorTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_WaterTemperature','Units',Units.Average_WaterTemperature);
ncwriteatt([data_path deployment(1:end-1) '_AD2CP.nc'],'Average_Range','Units',Units.Average_Range);



