function [a, phi] = find_principal_rotation(C21)
% find_principal_rotation
% From a DCM C21, return (a, phi) for the principal rotation (axis-angle).
%
% Uses:
%   cos(phi) = (tr(C) - 1)/2
%   a = (1/(2 sin(phi))) * [C32 - C23; C13 - C31; C21 - C12]
%
% Supports numeric or symbolic input.

    if isa(C21, 'sym')
        trC = trace(C21);
        phi = acos((trC - 1)/2);
        denom = 2*sin(phi);

        a = [C21(3,2) - C21(2,3);
             C21(1,3) - C21(3,1);
             C21(2,1) - C21(1,2)] / denom;

        n = norm(a);
        if ~isequal(n, sym(0))
            a = a / n;
        end
    else
        C = reshape(double(C21), 3, 3);
        trC = trace(C);
        phi = acos((trC - 1.0)/2.0);
        denom = 2*sin(phi);

        a = [C(3,2) - C(2,3);
             C(1,3) - C(3,1);
             C(2,1) - C(1,2)] / denom;

        n = norm(a);
        if n ~= 0
            a = a / n;
        end
    end
end