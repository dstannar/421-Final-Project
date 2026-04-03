clc
clear all

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