function [z] = cell_vert(pitch,roll,rng_cells)

% 
% % calculate a vertical displacment above the instrument for each adcp bin
% % adjusting for pitch and roll (in degrees)

% rng_cells is the vector of cell ranges found in the adcp .hdr file

junk  = rng_cells .* sin(deg2rad(90-pitch));
z = junk .* sin(deg2rad(90-roll));
clear junk

% figure
% plot(zeros(length(rng_cells),1),rng_cells,'bo')
% hold on
% grid
% plot(zeros(length(rng_cells),1),z,'r+')


return
