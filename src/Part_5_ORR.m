% Hailey Orr 421 Final Project Part 5
% control system design
% modify Part 2 

clc
clear all
close all

% givens in detumble ops 

w_b_ECI_0 = [-0.05
    0.03
    0.2];

d.mass = 640;
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
J = [d.Ixx, d.Ixy, d.Ixz; ...
       d.Ixy, d.Iyy, d.Iyz; ...
       d.Ixz, d.Iyz, d.Izz];

m = 640;

% given properties

mu = 398600;
h = 53335.2;
ecc = 0;
OMEGA = 0 * pi/180;
inc = deg2rad(98.43);
argofp = 0;
ta = 0;

a = (h^2 / mu) * (1 / 1 - ecc^2);
period = 2 * pi * sqrt(a^3 / mu);

T_d = [0; 0; 0]; % change so just Td is zero, but adding a commanded torque

% get position and velocity
[r_eci0, v_eci0] = COES2RV(h, inc, OMEGA, ecc, argofp, ta, mu);

% now need to get vectors into lvlh frame

v0 = norm(v_eci0);
r0 = norm(r_eci0);

% definititon of lvlh frame
z_lvlh = -r_eci0 / r0;
y_lvlh = - cross(r_eci0, v_eci0) / norm(cross(r_eci0, v_eci0));
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

% now add control law

Kp = zeros(3);
Kd = [0.2 0 0 
    0 0.2 0
    0 0 0.2];

% what is our commanded euler??
E_c = [0;0;0];
q_c = [E_c; sqrt(1 - norm(E_c)^2)];


tspan = period * 5; % wants 5 orbits

out = sim('ORR_Part5');

t = out.tout;
w = squeeze(out.w.signals.values);
E = squeeze(out.E.signals.values);
q = squeeze(out.q.signals.values);

% this one was being weird
M_c_sig = out.logsout.get('M_c');
M_c = squeeze(M_c_sig.Values.Data);



figure()
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

subplot(3, 1, 3)
plot(t, q)
xlabel('Time (secs)')
ylabel('q')
legend('q_1','q_2','q_3','q_4')
title('Quaternions')

figure()
plot(t, M_c)
xlabel('time (seconds')
ylabel('Torque (Nm)')
legend('T_x', 'T_y', 'T_z')
title('Thruster Torque for Detumble Phase')
% ------------------------
% ------------------------------------------------------------

% functions 

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

function quatconj = quatConjugate(q)

    quatconj = [-q(1:3);
        q(4)];

end

function w = quatMult(p, q)

    e_p = p(1:3);
    e_q = q(1:3);
    eta_p = p(4);
    eta_q = q(4);

    w = [eta_p*e_q + eta_q*e_p + vcross(e_p)*e_q;
        eta_p*eta_q - e_p'*e_q];
end