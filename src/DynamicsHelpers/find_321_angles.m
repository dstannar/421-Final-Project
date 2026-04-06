function [phi, theta, psi] = find_321_angles(C, degs)
% find_321_angles
% Extract 3-2-1 Euler angles from a DCM

    if nargin < 2
        degs = false;
    end

    phi = atan2(C(2,3), C(3,3));
    theta = -asin(C(1,3));
    psi = atan2(C(1,2), C(1,1));

    if degs
        phi = rad2deg(phi);
        theta = rad2deg(theta);
        psi = rad2deg(psi);
    end
end