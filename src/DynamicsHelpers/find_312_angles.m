function [psi, theta, phi] = find_312_angles(C, degs)
% find_312_angles
% Extract 3-1-2 Euler angles from a DCM
%
% Based on:
%   R = Cz(psi) * Cx(theta) * Cy(phi)

    if nargin < 2
        degs = false;
    end

    theta = asin(C(3,2));
    psi   = atan2(-C(1,2), C(2,2));
    phi   = atan2(-C(3,1), C(3,3));

    if degs
        psi = rad2deg(psi);
        theta = rad2deg(theta);
        phi = rad2deg(phi);
    end
end