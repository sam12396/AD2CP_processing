function M = CalcTiltMatrix(pitch, roll)

sinpp = sind(pitch);
cospp = cosd(pitch);
sinrr = sind(roll);
cosrr = cosd(roll);

M = [cospp -sinpp*sinrr  -cosrr*sinpp;...
       0       cosrr         -sinrr;  ...
     sinpp  sinrr*cospp  cospp*cosrr ]; 

end
