function C = DCM_from_principal(a, phi)
% DCM_from_principal
% Implements:
%   C = cos(phi) I + (1 - cos(phi)) a a^T - sin(phi) a_x
%
% where a is a unit vector and a_x = skew(a).
%
% Supports numeric or symbolic input.

    if isa(a, 'sym') || isa(phi, 'sym')
        a = a(:);
        n = norm(a);
        if ~isequal(n, sym(0))
            a = a / n;
        end
        I = sym(eye(3));
        C = cos(phi)*I + (1 - cos(phi))*(a*a.') - sin(phi)*skew(a);
    else
        a = a(:);
        n = norm(a);
        if n ~= 0
            a = a / n;
        end
        I = eye(3);
        aaT = a * a.';
        C = cos(phi)*I + (1 - cos(phi))*aaT - sin(phi)*skew(a);
    end
end