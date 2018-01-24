function e_plane = calcPlaneVectors(normal)
% calculates the plane vectors from a normal based on the criterion that
% the first plane vector is the normal, the second is the vector
% perpependicular to the normal and parallel to the XY-plane and the third
% vector is perpendicular to both
% EHarry October 2011, taken from a subfunction of makiFitPlane.m by
% KJaqaman

e_plane = zeros(3);
e_plane(:,1) = normal;
e_plane(:,2) = [-normal(2),normal(1),0]./sqrt(sum(normal(1:2).^2));
e_plane(:,3) = cross(e_plane(:,1),e_plane(:,2));

end