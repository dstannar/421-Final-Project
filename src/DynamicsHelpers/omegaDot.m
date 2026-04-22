function wdot = omegaDot(t, w, I, T, psinew)
% omegaDot
% Get time derivative of angular velocity

    wdot = I \ (T(:) - cross(w, I*w));
end