function dydt = quatDot(t, y, psinew, I, T)
% quatDot
% Get time derivative of quaternion and angular velocity

    w   = y(1:3);
    eps = y(4:6);
    eta = y(7);

    wdot   = omegaDot(t, w, I, T, psinew);
    epsdot = 0.5 * ((eta*eye(3) + skew(eps)) * w);
    etadot = -0.5 * (eps.' * w);

    dydt = [wdot; epsdot; etadot];
end