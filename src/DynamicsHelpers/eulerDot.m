function dydt = eulerDot(t, y, psinew, I, T)
% eulerDot
% Get time derivative of angular velocity and Euler angles

    phi = y(4);
    theta = y(5);
    psi = y(6);

    w = y(1:3);

    B = [1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
         0, cos(phi),            -sin(phi);
         0, sin(phi)/cos(theta),  cos(phi)/cos(theta)];

    angleDots = B * w;
    wdot = omegaDot(t, w, I, T, psinew);

    dydt = [wdot; angleDots];
end