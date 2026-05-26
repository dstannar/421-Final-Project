clear; close all; clc;

% get mass props (part 1)
%% DETUMBLE MODE (CUBE)
% mass, xyz dimensions, center of mass, inertia matrix

d.mass = 640; % [kg]
d.x = 2; % [m]
d.y = 2; % [m]
d.z = 2; % [m]
d.COM = [0, 0, 0]; 

d.Ixx = (1/12)*d.mass*(d.y^2 + d.z^2);
d.Iyy = (1/12)*d.mass*(d.x^2 + d.z^2);
d.Izz = (1/12)*d.mass*(d.x^2 + d.y^2);
d.Ixy = 0;
d.Iyz = 0;
d.Ixz = 0;

d.I = [d.Ixx, d.Ixy, d.Ixz; ...
       d.Ixy, d.Iyy, d.Iyz; ...
       d.Ixz, d.Iyz, d.Izz];

%% NORMAL MODE (CUBE + SENSOR + 2X PANEL)
% mass, xyz dimensions, center of mass, inertia matrix

% bus
n.bus.mass = 500; % [kg]
n.bus.x = 2; % [m]
n.bus.y = 2; % [m]
n.bus.z = 2; % [m]
n.bus.COM = [0, 0, 0];

n.bus.Ixx = (1/12)*n.bus.mass*(n.bus.y^2 + n.bus.z^2);
n.bus.Iyy = (1/12)*n.bus.mass*(n.bus.x^2 + n.bus.z^2);
n.bus.Izz = (1/12)*n.bus.mass*(n.bus.x^2 + n.bus.y^2);
n.bus.Ixy = 0;
n.bus.Iyz = 0;
n.bus.Ixz = 0;
n.bus.I = [n.bus.Ixx, n.bus.Ixy, n.bus.Ixz; ...
           n.bus.Ixy, n.bus.Iyy, n.bus.Iyz; ...
           n.bus.Ixz, n.bus.Iyz, n.bus.Izz];

% sensor
n.sens.mass = 100; % [kg]
n.sens.x = 0.25; % [m]
n.sens.y = 0.25; % [m]
n.sens.z = 1; % [m]
n.sens.COM = [0, 0, 1.5];

n.sens.Ixx = (1/12)*n.sens.mass*(n.sens.y^2 + n.sens.z^2);
n.sens.Iyy = (1/12)*n.sens.mass*(n.sens.x^2 + n.sens.z^2);
n.sens.Izz = (1/12)*n.sens.mass*(n.sens.x^2 + n.sens.y^2);
n.sens.Ixy = 0;
n.sens.Iyz = 0;
n.sens.Ixz = 0;
n.sens.I = [n.sens.Ixx, n.sens.Ixy, n.sens.Ixz; ...
            n.sens.Ixy, n.sens.Iyy, n.sens.Iyz; ...
            n.sens.Ixz, n.sens.Iyz, n.sens.Izz];

% left panel (-y)
n.panL.mass = 20; % [kg]
n.panL.x = 2; % [m]
n.panL.y = 3; % [m]
n.panL.z = 0.05; % [m]
n.panL.COM = [0, -2.5, 0];

n.panL.Ixx = (1/12)*n.panL.mass*(n.panL.y^2 + n.panL.z^2);
n.panL.Iyy = (1/12)*n.panL.mass*(n.panL.x^2 + n.panL.z^2);
n.panL.Izz = (1/12)*n.panL.mass*(n.panL.x^2 + n.panL.y^2);
n.panL.Ixy = 0;
n.panL.Iyz = 0;
n.panL.Ixz = 0;
n.panL.I = [n.panL.Ixx, n.panL.Ixy, n.panL.Ixz; ...
            n.panL.Ixy, n.panL.Iyy, n.panL.Iyz; ...
            n.panL.Ixz, n.panL.Iyz, n.panL.Izz];

% right panel (+y)
n.panR.mass = 20; % [kg]
n.panR.x = 2; % [m]
n.panR.y = 3; % [m]
n.panR.z = 0.05; % [m]
n.panR.COM = [0, 2.5, 0];

n.panR.Ixx = (1/12)*n.panR.mass*(n.panR.y^2 + n.panR.z^2);
n.panR.Iyy = (1/12)*n.panR.mass*(n.panR.x^2 + n.panR.z^2);
n.panR.Izz = (1/12)*n.panR.mass*(n.panR.x^2 + n.panR.y^2);
n.panR.Ixy = 0;
n.panR.Iyz = 0;
n.panR.Ixz = 0;
n.panR.I = [n.panR.Ixx, n.panR.Ixy, n.panR.Ixz; ...
            n.panR.Ixy, n.panR.Iyy, n.panR.Iyz; ...
            n.panR.Ixz, n.panR.Iyz, n.panR.Izz];

% total satellite
n.sat.mass = n.bus.mass + n.sens.mass + n.panL.mass + n.panR.mass; % [kg]
n.sat.COM = ((n.bus.mass*n.bus.COM) + (n.sens.mass*n.sens.COM) + (n.panL.mass*n.panL.COM) + (n.panR.mass*n.panR.COM))/n.sat.mass;

% sum all inertia matrices with some parallel axis theorem
n.sat.I = (n.bus.I  - n.bus.mass *(vectX(n.sat.COM - n.bus.COM))*(vectX(n.sat.COM - n.bus.COM))) + ...
          (n.sens.I - n.sens.mass*(vectX(n.sat.COM - n.sens.COM))*(vectX(n.sat.COM - n.sens.COM))) + ...
          (n.panL.I - n.panL.mass*(vectX(n.sat.COM - n.panL.COM))*(vectX(n.sat.COM - n.panL.COM))) + ...
          (n.panR.I - n.panR.mass*(vectX(n.sat.COM - n.panR.COM))*(vectX(n.sat.COM - n.panR.COM)));

% cross operator function necessary for parallel axis theorem
function vCross = vectX(v)
    vCross = [ 0,    -v(3),  v(2)
              v(3),    0,   -v(1)
             -v(2),   v(1),   0];
end

%% RESPONSES TO FINAL PROJECT PART 1
fprintf("Detumble total mass = %.0f kg\n", d.mass)
fprintf("Detumble COM =\n")
disp(d.COM)
fprintf("Detumble inertia matrix =\n")
disp(d.I)

fprintf("Normal total mass = %.0f kg\n", n.sat.mass)
fprintf("Normal COM =\n")
disp(n.sat.COM)
fprintf("Normal inertia matrix =\n")
disp(n.sat.I)

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

%% Part 5 - Detumble

w_0 = w0_b_eci;

% Detumble control law:
% Mc = -Kd*w_b_ECI

% Simulate for five orbits
tspan = 5*orbital_period;

J = n.sat.I;

% only detumble, don't care what attitude but need to not error 
Kp = zeros(3);
q_C = [0; 0; 0; 1];

%% Set Kd
Mp_reg = .1;
ts = 10^4;

zeta = sqrt(log(Mp_reg)^2/(pi^2 + log(Mp_reg)^2));

wn = log(0.02*sqrt(1-zeta^2))/-zeta/ts;

beta = atan(sqrt(1-zeta^2)/zeta);
tr = (pi-beta)/wn/sqrt(1-zeta^2);

Kd = J/ts

% not modeling disturbance torque
T_d = [0; 0; 0];

% Run Simulink
out = sim('SS_Part5');

q_b_ECI = out.get('q_b_ECI');
w_b_ECI = out.get('w_b_ECI');
M_c     = out.get('M_c');

%% Extract Simulink outputs

% q_b_ECI, w_b_ECI, and M_c
t_q = q_b_ECI.time;
q_vals = squeeze(q_b_ECI.signals.values);

t_w = w_b_ECI.time;
w_vals = squeeze(w_b_ECI.signals.values);

t_M = M_c.time;
M_vals = squeeze(M_c.signals.values);

%% Convert quaternion history to 3-2-1 Euler angles

Eul_b_ECI = zeros(length(t_q),3);
q_vals = squeeze(q_b_ECI.signals.values);

for k = 1:length(t_q)
    C_b_ECI = quat_to_DCM(q_vals(k,:).');
    [phi, theta, psi] = find_321_angles(C_b_ECI, false);
    Eul_b_ECI(k,:) = [phi, theta, psi];
end

%% Plots

Eul_b_ECI = unwrap(Eul_b_ECI);

figure
subplot(3,1,1)
plot(t_q, rad2deg(Eul_b_ECI), 'LineWidth', 1.2)
xlabel('Time (secs)')
ylabel('Euler Angles (deg)')
legend('\phi','\theta','\psi')
title('3-2-1 Euler Angles: Body to ECI')
grid on

subplot(3,1,2)
plot(t_q, q_vals, 'LineWidth', 1.2)
xlabel('Time (secs)')
ylabel('Quaternion')
legend('\epsilon_x','\epsilon_y','\epsilon_z','\eta')
title('Quaternion: Body to ECI')
grid on

subplot(3,1,3)
plot(t_w, w_vals, 'LineWidth', 1.2)
xlabel('Time (secs)')
ylabel('Body Rates (rad/s)')
legend('\omega_x','\omega_y','\omega_z')
title('Angular Velocity in Body Components')
grid on

figure
plot(t_M, M_vals, 'LineWidth', 1.2)
xlabel('Time (secs)')
ylabel('Torque (N-m)')
legend('M_x','M_y','M_z')
title('Commanded Torque in Body Frame')
grid on