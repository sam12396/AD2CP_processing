function heading = CalcMidlifeHeading(HxHyHz_sensor, pitch, roll, orientation)

% *** Start of magnetometer recalibration ***


% *** End of magnetometer recalibration ***

% Heading calculation
M_tilt = CalcTiltMatrix(pitch, roll);

if orientation == 0
    % Upwards looking instrument
    HxHyHz_inst = HxHyHz_sensor;
else
    % Downwards looking instrument (rotation around x-axis)
    HxHyHz_inst(1)   =  HxHyHz_sensor(1);
    HxHyHz_inst(2:3) = -HxHyHz_sensor(2:3);
end

% Calculate the magnetic vector in the Earth aligned coordinate system
HxHyHz_Earth = M_tilt*HxHyHz_inst';

% Calculate the heading
heading = atan2(HxHyHz_Earth(2),HxHyHz_Earth(1))*180/pi;
if heading < 0
    heading = heading + 360;
end

end