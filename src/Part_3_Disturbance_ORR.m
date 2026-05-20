% Hailey Orr AERO 421, Final Project Part 3

clc
clear all 

%% Part 1 ---------------------------------------------------------------

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



%% RESPONSES TO FINAL PROJECT PART 1
fprintf("Normal total mass = %.0f kg\n", n.sat.mass)
fprintf("Normal COM =\n")
disp(n.sat.COM)
fprintf("Normal inertia matrix =\n")
disp(n.sat.I)

%% Part 2 --------------------------------------------------------------------------
%% UPDATE
% satellite bus 
Areas = 4*ones(6,1);
normals = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
cps = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

% Append geometric properties for Solar Panel 1
% Append geometric properties for Solar Panel 2
% Append geometric properties for Sensor
% now subtract the center of mass to get the location of the rho vectors
% with respect to the center of mass
% Now build the matrix

surfaceProperties = [Areas cps normals];

%% Part 3 -----------------------------------------------------------------

% magnetic field model
JD_0 = 2460390;

% Spacecraft Orbit Properties
mu = 398600; % km^3/s^2
h = 53335.2; % km^2/s
e = 0; % none
Omega = 0*pi/180; % radians
inclination = 98.43*pi/180; % radians
omega = 0*pi/180; % radians
nu = 0*pi/180; % radians
a = h^2/mu/(1 - e^2);
orbital_period = 2*pi*sqrt(a^3/mu);

% Set/Compute initial conditions
% intial orbital position and velocity
[r_ECI_0, v_ECI_0]= COES2RV(h, inclination, Omega, e, omega, nu, mu);

% No external command Torque
T_c = [0; 0; 0]; % Nm

% now need to get vectors into lvlh frame

v0 = norm(v_ECI_0);
r0 = norm(r_ECI_0);

% definititon of lvlh frame
z_lvlh = -r_ECI_0 / r0;
y_lvlh = - cross(r_ECI_0, v_ECI_0) / norm(cross(r_ECI_0, v_ECI_0));
x_lvlh = cross(y_lvlh, z_lvlh);

% initial euler
phi0 = 0;
theta0 = 0;
psi0 = 0;
E_b_lvlh_0 = [phi0; theta0; psi0];

% quarternions
q_b_lvlh_0 = [0; 0; 0; 1];

% compute initial rotation matrices
C_lvlh_eci_0 = [x_lvlh'; y_lvlh'; z_lvlh'];

C_b_lvlh_0 = [1 0 0
    0 1 0
    0 0 1];

C_b_eci_0 = C_lvlh_eci_0 * C_b_lvlh_0;

E_b_ECI_0 = C2euler(C_b_eci_0);

q_b_ECI_0 = C2quat(C_b_eci_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

%% Part 4 - Simulate Results

n_revs = 5; %revs
tspan = n_revs * orbital_period;
out = sim('FP_Solutions_Part_3_disturbance');


%% Part 5 - Plot Results
% Plot Angular Velocities, Euler Angles and Quaternions
% Plot Disturbance torques in F_b

t = out.tout;
w = squeeze(out.w.signals.values);
E = squeeze(out.E.signals.values);
q = squeeze(out.q.signals.values);

figure
subplot(3,1,1)
plot(t, w)
ylabel('\omega (rad/s)')
xlabel('Time (secs)')
legend('\omega_x','\omega_y','\omega_z')
grid on
title('Angular Velocities')

subplot(3,1,2)
plot(t, E)
ylabel('Euler Angles (radians)')
xlabel('Time (secs)')
legend('\phi','\theta','\psi')
grid on
title('Euler Angles')

subplot(3,1,3)
plot(t, q)
xlabel('Time (secs)')
ylabel('q')
legend('q_1','q_2','q_3','q_4')
title('Quaternions')



%% Functions 

% cross operator function necessary for parallel axis theorem
function vCross = vectX(v)
    vCross = [ 0,    -v(3),  v(2)
              v(3),    0,   -v(1)
             -v(2),   v(1),   0];
end

function [r, v] = COES2RV(h, inc, RAAN, ecc, omega, TA, mu)

% make sure everything is in radians!!

    rp = (h^2/mu) * (1/(1 + ecc * cos(TA))) * (cos(TA) * [1;0;0] + sin(TA) * [0;1;0]);
    vp = (mu/h) * (-sin(TA) * [1;0;0] + (ecc + cos(TA)) * [0;1;0]);

    R3_RAAN = [cos(RAAN) sin(RAAN) 0
         -sin(RAAN) cos(RAAN) 0
         0 0 1];

    R1_inc = [1 0 0
     0 cos(inc) sin(inc)
     0 -sin(inc) cos(inc)];

    R3_omega = [cos(omega) sin(omega) 0
     -sin(omega) cos(omega) 0
     0 0 1];

    Q_pX = (R3_omega * R1_inc * R3_RAAN)';

    r = Q_pX * rp;
    v = Q_pX * vp;

end



function euler = C2euler(C21)

    phi = atan(C21(2,3) / C21(3,3));
    theta = -asin(C21(1,3));
    psi = atan(C21(1,2) / C21(1,1));

    euler = [phi; theta; psi];

end


function quaternion = C2quat(C21)
   
    phi = acos((trace(C21) - 1) / 2);
    eta0 = (trace(C21) + 1)^(0.5) / 2;

    a10 = (C21(2,3) - C21(3,2)) / (4 * eta0);
    a20 = (C21(3,1) - C21(1,3)) / (4 * eta0);
    a30 = (C21(1,2) - C21(2,1)) / (4 * eta0);

    o = sin(phi / 2); % random variable

    epsilon0 = [a10 * o; a20 * o; a30 * o]; 
    quaternion = [epsilon0; eta0];

end