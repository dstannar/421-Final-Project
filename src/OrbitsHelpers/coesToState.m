function [rvec, vvec] = coesToState(hmag, ecc, TA, raan, inc, argp, mu)
% coesToState
% Convert classical orbital elements to inertial state vectors
%
% Inputs:
%   hmag : specific angular momentum magnitude [km^2/s]
%   ecc  : eccentricity
%   TA   : true anomaly [rad]
%   raan : right ascension of ascending node [rad]
%   inc  : inclination [rad]
%   argp : argument of perigee [rad]
%   mu   : gravitational parameter [km^3/s^2]
%
% Outputs:
%   rvec : 3x1 position vector in inertial frame [km]
%   vvec : 3x1 velocity vector in inertial frame [km/s]

    p = hmag^2 / mu;
    rmag = p / (1 + ecc*cos(TA));

    % State in PQW
    r0 = rmag * [cos(TA); sin(TA); 0.0];
    v0 = (mu/hmag) * [-sin(TA); ecc + cos(TA); 0.0];

    % Rotation matrices
    R3_W = [ cos(raan),  sin(raan), 0;
            -sin(raan),  cos(raan), 0;
                     0,          0, 1];

    R1_i = [1,        0,         0;
            0, cos(inc),  sin(inc);
            0,-sin(inc),  cos(inc)];

    R3_w = [ cos(argp),  sin(argp), 0;
            -sin(argp),  cos(argp), 0;
                     0,          0, 1];

    Q_peri_ECI = R3_w * R1_i * R3_W;
    Q_ECI_peri = Q_peri_ECI.';

    rvec = Q_ECI_peri * r0;
    vvec = Q_ECI_peri * v0;
end