function params = getMoreOrbitalParams(hmag, ecc, TA, raan, inc, argp, mu)
% getMoreOrbitalParams
% Compute a broad set of useful orbital parameters from COEs

% Inputs:
%   hmag : specific angular momentum magnitude [km^2/s]
%   ecc  : eccentricity
%   TA   : true anomaly [rad]
%   raan : RAAN [rad]
%   inc  : inclination [rad]
%   argp : argument of perigee [rad]
%   mu   : gravitational parameter [km^3/s^2]
%
% Output:
%   params : struct with many derived orbital quantities

    % Semilatus rectum
    p = hmag^2 / mu;

    % Radius at current true anomaly
    rmag = p / (1 + ecc*cos(TA));

    % Periapsis / apoapsis radii
    r_per = hmag^2 / mu * (1 / (1 + ecc*cos(0)));
    r_apo = hmag^2 / mu * (1 / (1 + ecc*cos(pi)));

    % Semimajor axis
    sma = 0.5 * (r_per + r_apo);

    % Period (elliptic case, same as Python)
    period = 2*pi / sqrt(mu) * sma^(3/2);

    % State vectors from COEs
    [rvec, vvec] = coesToState(hmag, ecc, TA, raan, inc, argp, mu);

    % Magnitudes and radial/transverse velocity
    vmag = norm(vvec);
    v_r = dot(rvec, vvec) / rmag;
    v_t = hmag / rmag;

    % Flight path angle
    FPA = atan2(v_r, v_t);

    % Specific orbital energy
    energy = vmag^2 / 2 - mu / rmag;

    % Eccentric anomaly and mean anomaly
    EA = 2 * atan2( sqrt(1 - ecc) * sin(TA/2), ...
                    sqrt(1 + ecc) * cos(TA/2) );
    EA = mod(EA, 2*pi);

    MA = mod(EA - ecc*sin(EA), 2*pi);

    % Mean motion
    n = sqrt(mu / sma^3);

    % Orbit geometry
    b = sma * sqrt(1 - ecc^2);      % semiminor axis
    c = sma * ecc;                  % focal distance

    % Perifocal vectors
    evec_pqw = [ecc; 0; 0];
    hvec_pqw = [0; 0; hmag];

    % Inertial-frame h and e vectors
    hvec = cross(rvec, vvec);
    eccvec = (1/mu) * (((vmag^2 - mu/rmag) * rvec) - (rmag * v_r * vvec));

    % Speeds at apsides
    v_per = mu / hmag * (1 + ecc);
    v_apo = mu / hmag * (1 - ecc);

    % Orbital plane normal unit vector
    hhat = hvec / norm(hvec);

    params = struct();
    params.hmag   = hmag;
    params.ecc    = ecc;
    params.TA     = TA;
    params.raan   = raan;
    params.inc    = inc;
    params.argp   = argp;

    params.p      = p;
    params.sma    = sma;
    params.b      = b;
    params.c      = c;
    params.rmag   = rmag;
    params.r_per  = r_per;
    params.r_apo  = r_apo;

    params.period = period;
    params.n      = n;
    params.energy = energy;

    params.EA     = EA;
    params.MA     = MA;

    params.rvec   = rvec;
    params.vvec   = vvec;
    params.vmag   = vmag;
    params.v_r    = v_r;
    params.v_t    = v_t;
    params.FPA    = FPA;

    params.hvec   = hvec;
    params.hhat   = hhat;
    params.eccvec = eccvec;
    params.evec_pqw = evec_pqw;
    params.hvec_pqw = hvec_pqw;

    params.v_per  = v_per;
    params.v_apo  = v_apo;
end