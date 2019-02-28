function [O_ls,G_ls,bnew,C] = inversion_leastSquare_sparse_2019(U,V,Z,dz,uv_daverage)
%[O_ls,G_ls,bnew,C] = inversion_leastSquare_sparse_2019(U,V,Z,dz,uv_daverage)
% Prior to running this remove all bad data points above, below, or in the
% glider dive.
% dz is desired vertical resolution, but should not be much smaller than
% the bin length
% U is measured east-west velocities from ADCP
% V is measured north-south velocities from ADCP
% Z is the measurement depths of U and V
% uv_daverage is depth averaged velocity
% O_ls is the ocean velocity profile
% G_ls is the glider velocity profile
% bnew are the bin centers for the point in the profiles
% C is the constant used in the constraint equation

% This function's purpose is to take velocity measurements from a glider
% mounted ADCP and separate the ocean velocity and glider velocity from
% each of these measurements.
%% To-do list
%  1) Chop off any measurements from below max glider depth or pressure
%       This might happen before the data gets here but I have to check
%  2) Need to add curvature-minimizing smoothness constraint

%% QA/QC Step 1 Check that all depths have data associated with them

jj= isnan(U+V+Z);
U(jj)=NaN;
V(jj)=NaN;
Z(jj)=NaN;

clear jj

%% Visbeck 2002 variables

% vars from Visbeck (2002) page 800
% eqn(19) Knowns
% Maximum number of observations (nd) is given by the number of velocity
% estimates per ping (nbin) times the number of profiles per cast (nt)
nbin = size(U,2);     % number of programmed ADCP bins per individual profile
nt = size(U,1);       % number of individual velocity profiles
nd = nbin*nt;         % G dimension (1) 


% Check to see if orientation of UVZ is correct
if nbin>nt
    disp('nbin>nt Check dimensions of the inputs are [profiles bins]')
end

% eqn(21) Unknowns are the sum of the ocean velocities plus the CTD
% velocities for each ensemble
H = floor(max(Z(:))); % max measurement depth for this particular file 
bedge = 0:dz:H;       % define the edges of the bins

%% Check that each bin has data in it

meas=zeros(length(bedge)-1,1); %pre allocate memory
for k=1:length(bedge)-1
    % Creates index of Z values that fall inside the bin edges
    ii= Z >= bedge(k) & Z < bedge(k+1);
    meas(k)=sum(sum(ii));
    clear ii;
end

bnew = bedge(1:end-1) + dz/2; % These are bin centers

if size(bnew)==[1,1]
    error('Desired profile resolution, dz, is too large for depth of profile')
end


% Chop off the top of the profile that does not have data 
ii=meas > 0;                %Finds all bins with data
ij = find(ii, 1, 'first');  %Finds first bin with data
bnew(1:ij-1)=[];            %Removes all bins above the first with data
z1=bnew(1);                 %z1 is the depth of the midpoint of the first bin with data

%% Create and populate G

nz=length(bnew); % number of ocean velocities desired
nm = nz + nt; % G dimension (2), number of unknowns

% Let's build the corresponding coefficient matrix G 
G = zeros(nd,nm);

for ii = 1:nt % number of adcp profiles per cast
    for jj = 1:nbin % number of measured bins per profile
        
  % We will skip this (leave it = 0) if U(nt(ii),nbin(jj)) = nan
        if isfinite(U(ii,jj))
            
            % This section will take care of the Uctd part of the matrix
            G((nbin*(ii-1))+jj,ii) = 1;        
        
            % This will fill in the Uocean part
            % We loop through all Z members and place them in the proper 
            % G matrix location
            
            % Find the difference between all bin centers and the current Z
            % value
            dx=abs(bnew-Z(ii,jj));      
            % Find the minimum of these differences
            minx = min(dx);
            % Finds bnew index of the first match of Z and bnew
            idx = find((dx-minx)==0,1);  %% Need to get rid of 'find'      
            G((nbin*(ii-1))+jj,nt+idx) = 1;       
            clear dx minx idx
        end % of isfinite(U) switch
        
    end % of jj
end % of ii
clear ii jj

%% Reshape U and V into the d format

d_u = U.'; % transpose
d_u = d_u(:); % and form a column
    
d_v = V.'; 
d_v = d_v(:);

%% Need to calculate C (Todd et al. 2017) based on our inputs 
% This creates a row that has the same # of columns as G. The elements
% of the row follow the trapezoid rule which is used because of the
% extension of the first bin with data to the surface. The last entry of
% the row corresponds to the max depth reached by the glider, any bins
% below that should have already been removed.
constraint=[zeros(1,nt),  z1/2 (z1/2+dz/2) repmat(dz,1,nz-3) dz/2];

% To find see we use the equation of the norm and set norm=1 because we
% desire unity. The equation requires we take the sum of the squares of the
% entries in constraint.
sqr_constraint=constraint.*constraint;
sum_sqr_constraint=sum(sqr_constraint(:));

%Then we can solve for the value of C needed to maintain unity 
C=H*(1/sqrt(sum_sqr_constraint));


%% Applying the depth averaged velocity constraint from Todd et al., (2011/2017)
% Add term to the d arrays
du = [d_u; C*uv_daverage(1)]; 
dv = [d_v; C*uv_daverage(2)];

d = complex(du,dv);
Gstar = [G; (C/H)*constraint];

%% INVERSION OF THE SPARSE MATRIX
 
Gs=sparse(Gstar);
ms = ((Gs'*Gs)\Gs')*d; % Least square inverse solution


O_ls = ms(nt+1:end); % Ocean velocity
G_ls = ms(1:nt); % Glider velocity
clear Gs ms

return