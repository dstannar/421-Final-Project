% Hailey Orr 421 Final Project Part 2
% torque free motion

clc
clear all

% givens in normal ops

w_b_ECI_0 = [0.001
    -0.001
    0.002];

J = [812.0396         0         0
         0  545.3729         0
         0         0  627.7083];

m = 640;

com = [0 0 0.2344];

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

T = [0; 0; 0];

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

q_b_ECI_0 = C2quart(C_b_eci_0);

tspan = period;

out = sim('ORR_Part2');




% ------------------------------------------------------------------------------------

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


function quarternion = C2quart(C21)
   
    phi = acos((trace(C21) - 1) / 2);
    eta0 = (trace(C21) + 1)^(0.5) / 2;

    a10 = (C21(2,3) - C21(3,2)) / (4 * eta0);
    a20 = (C21(3,1) - C21(1,3)) / (4 * eta0);
    a30 = (C21(1,2) - C21(2,1)) / (4 * eta0);

    o = sin(phi / 2); % random variable

    epsilon0 = [a10 * o; a20 * o; a30 * o]; 
    quarternion = [epsilon0; eta0];

end