function coes = stateToCOEs(rvec, vvec, mu)
% stateToCOEs
% Convert inertial state vectors to classical orbital elements (COEs)
%
% Inputs:
%   rvec : 3x1 or 1x3 position vector [km]
%   vvec : 3x1 or 1x3 velocity vector [km/s]
%   mu   : gravitational parameter [km^3/s^2]
%
% Output:
%   coes : struct containing:
%       .hmag
%       .inc
%       .raan
%       .ecc
%       .argp
%       .TA
%       .r_per
%       .r_apo
%       .sma
%       .period
%       .energy
%       .EA
%       .MA
%       .hvec
%       .v_r
%       .v_t
%       .FPA
%       .vmag
%       .rmag
%       .eccvec

    rvec = rvec(:);
    vvec = vvec(:);
    k_hat = [0; 0; 1];

    rmag = norm(rvec);
    vmag = norm(vvec);
    v_r  = dot(rvec, vvec) / rmag;

    h = cross(rvec, vvec);
    hmag = norm(h);
    inc = acos(h(3) / hmag);

    Nvec = cross(k_hat, h);
    Nmag = norm(Nvec);

    % RAAN
    raan_raw = acos(Nvec(1) / Nmag);
    if Nvec(2) >= 0
        raan = raan_raw;
    else
        raan = 2*pi - raan_raw;
    end

    % Eccentricity vector
    evec = (1/mu) * ( ((vmag^2 - mu/rmag) * rvec) - (rmag * v_r * vvec) );
    ecc = norm(evec);

    % Argument of perigee
    arg_per_raw = acos(dot(Nvec, evec) / (Nmag * ecc));
    if evec(3) >= 0
        argp = arg_per_raw;
    else
        argp = 2*pi - arg_per_raw;
    end

    % True anomaly
    ta_raw = acos(dot(evec, rvec) / (ecc * rmag));
    if v_r >= 0
        TA = ta_raw;
    else
        TA = 2*pi - ta_raw;
    end

    % More orbital parameters (same math as Python)
    r_per = hmag^2 / mu * (1 / (1 + ecc*cos(0)));
    r_apo = hmag^2 / mu * (1 / (1 + ecc*cos(pi)));
    sma = 0.5 * (r_per + r_apo);

    period = 2*pi / sqrt(mu) * sma^(3/2);
    energy = vmag^2 / 2 - mu / rmag;

    EA = 2 * atan2( sqrt(1 - ecc) * sin(TA/2), ...
                    sqrt(1 + ecc) * cos(TA/2) );
    EA = mod(EA, 2*pi);

    MA = mod(EA - ecc*sin(EA), 2*pi);

    v_t = hmag / rmag;
    FPA = atan2(v_r, v_t);

    coes = struct();
    coes.hmag   = hmag;
    coes.inc    = inc;
    coes.raan   = raan;
    coes.ecc    = ecc;
    coes.argp   = argp;
    coes.TA     = TA;
    coes.r_per  = r_per;
    coes.r_apo  = r_apo;
    coes.sma    = sma;
    coes.period = period;
    coes.energy = energy;
    coes.EA     = EA;
    coes.MA     = MA;
    coes.hvec   = h;
    coes.v_r    = v_r;
    coes.v_t    = v_t;
    coes.FPA    = FPA;
    coes.vmag   = vmag;
    coes.rmag   = rmag;
    coes.eccvec = evec;
end