function [ psi, theta, phi ] = eulerAnglesFromRotMat( rotMat )
%EULERANGLESFROMROTMAT gets euler roation angles from a rotation matrix,
%algorithm from https://truesculpt.googlecode.com/hg-history/38000e9dfece971460473d5788c235fbbe82f31b/Doc/rotation_matrix_to_euler.pdf
%   EHarry Dec 2012

if (abs(rotMat(3,1)) ~= 1)
    theta = -asin(rotMat(3,1));    
    phi = atan2(rotMat(3,2)/cos(theta),rotMat(3,3)/cos(theta));
    psi = atan2(rotMat(2,1)/cos(theta),rotMat(1,1)/cos(theta));
else
    psi = 0;
    if (rotMat(3,1) == -1)
        theta = pi/2;
        phi = atan2(rotMat(1,2),rotMat(1,3));
    else
        theta = -pi/2;
        phi = atan2(-rotMat(1,2),-rotMat(1,3));
    end
end

end

