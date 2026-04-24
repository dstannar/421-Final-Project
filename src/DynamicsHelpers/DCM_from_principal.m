function C = DCM_from_principal(a, phi)
% DCM_from_principal
% Implements:
%   C = cos(phi) I + (1 - cos(phi)) a a^T - sin(phi) a_x
%
% where a is a unit vector and a_x = skew(a).
%

    a = a(:);
    n = norm(a);
    if n ~= 0
        a = a / n;
    end
    I = eye(3);
    aaT = a * a.';
    C = cos(phi)*I + (1 - cos(phi))*aaT - sin(phi)*skew(a);
end