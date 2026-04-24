clear; close all; clc;

% get mass props
Part1_MassProps;

% load constants
orbitsConstants;

% orbit props
hmag = 53335.2;          % km^2/s
ecc  = 0;                % unitless
raan = 0*pi/180;         % rad
inc  = 98.43*pi/180;     % rad
argp = 0*pi/180;         % rad
TA   = 0*pi/180;         % rad

% initial body ang velocity wrt ECI, in body components
w0_b_eci = [0.001; -0.001; 0.002];   % rad/s

% torque free
T = [0;0;0];

% initial state vector in ECI
[r0_eci, v0_eci] = coesToState(hmag, ecc, TA, raan, inc, argp, muE);

% orbit period
orbParams = getMoreOrbitalParams(hmag, ecc, TA, raan, inc, argp, muE);
orbital_period = orbParams.period;

% initial attitude aligned with LVLH
phi0   = 0;
theta0 = 0;
psi0   = 0;
Eul0_b_lvlh = [phi0; theta0; psi0];

% initial quaternion relating body to LVLH
q0_b_lvlh = [0;0;0;1];

% DCM relating LVLH to ECI
% rows are LVLH basis vectors written in ECI components
C_LVLH_ECI_0 = LVLH_from_ECI(r0_eci, v0_eci);

% DCM relating body to LVLH
C_b_LVLH_0 = build_eulerDCM(3, 2, 1, phi0, theta0, psi0, false);

% DCM relating body to ECI
C_b_ECI_0 = C_b_LVLH_0 * C_LVLH_ECI_0;

% initial Euler angles relating body to ECI using 3-2-1 convention
[phi_b_eci_0, theta_b_eci_0, psi_b_eci_0] = find_321_angles(C_b_ECI_0, false);
Eul0_b_eci = [phi_b_eci_0; theta_b_eci_0; psi_b_eci_0];

% initial quaternion relating body to ECI
q0_b_eci = DCM_to_quat(C_b_ECI_0);

% simulate for one orbit
tspan = [0 orbital_period];

% USE SIMULINK
out = sim('SS_Part2.slx');