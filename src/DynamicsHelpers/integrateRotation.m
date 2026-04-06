function [QuatIVP, EulerIVP] = integrateRotation(phi, theta, psi, eps, eta, w, tspan, Ip, Tc_p, psinew)
% integrateRotation
% Integrates rotational dynamics and kinematics using Euler angles and quaternions
%
% Outputs:
%   QuatIVP  : struct with fields .t and .y
%   EulerIVP : struct with fields .t and .y

    yAngles = [phi; theta; psi];
    y0E = [w(:); yAngles(:)];
    y0Q = [w(:); eps(:); eta];

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

    [tE, yE] = ode45(@(t,y) eulerDot(t, y, psinew, Ip, Tc_p), tspan, y0E, opts);
    [tQ, yQ] = ode45(@(t,y) quatDot(t, y, psinew, Ip, Tc_p),  tspan, y0Q, opts);

    EulerIVP.t = tE;
    EulerIVP.y = yE.';

    QuatIVP.t = tQ;
    QuatIVP.y = yQ.';
end