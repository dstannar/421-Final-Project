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
%out = sim('Final_Project_Part_2_Simulink');

% use integrateRotation function
[QuatIVP, EulerIVP] = integrateRotation( ...
    phi_b_eci_0, theta_b_eci_0, psi_b_eci_0, ...
    q0_b_eci(1:3), q0_b_eci(4), w0_b_eci, ...
    tspan, n.sat.I, T, 0);

% unpack results
tE = EulerIVP.t;
yE = EulerIVP.y.'; % columns: [wx wy wz phi theta psi]

tQ = QuatIVP.t;
yQ = QuatIVP.y.'; % columns: [wx wy wz e1 e2 e3 eta]

w_hist = yE(:,1:3);
eul_hist = yE(:,4:6);
quat_hist = yQ(:,4:7);

% plots

% Euler angles
figure;
plot(tE, eul_hist(:,1), 'DisplayName', '\phi'); hold on;
plot(tE, eul_hist(:,2), 'DisplayName', '\theta');
plot(tE, eul_hist(:,3), 'DisplayName', '\psi');
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Euler Angles: Body to ECI');
legend();
grid on;

% Quaternion
figure;
plot(tQ, quat_hist(:,1), 'DisplayName', 'e_1'); hold on;
plot(tQ, quat_hist(:,2), 'DisplayName', 'e_2');
plot(tQ, quat_hist(:,3), 'DisplayName', 'e_3');
plot(tQ, quat_hist(:,4), 'DisplayName', '\eta');
xlabel('Time (s)');
ylabel('Quaternion Component');
title('Quaternion: Body to ECI');
legend();
grid on;

% Angular velocity
figure;
plot(tE, w_hist(:,1), 'DisplayName', '\omega_1'); hold on;
plot(tE, w_hist(:,2), 'DisplayName', '\omega_2');
plot(tE, w_hist(:,3), 'DisplayName', '\omega_3');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
title('Body Angular Velocity');
legend();
grid on;