function C = quat_to_DCM(q)
% quat_to_DCM
% Quaternion to DCM using the convention q = [epsilon; eta],
% where epsilon = [e1; e2; e3] is the vector part and eta is the scalar.
%
%   C = (2*eta^2 - 1) I + 2*epsilon*epsilon^T - 2*eta*(epsilon_x)
%
% Requires unit quaternion; this function normalizes q.
% Supports numeric or symbolic input.

    if isa(q, 'sym')
        q = q(:);

        n = sqrt(q.'*q);
        if ~isequal(n, sym(0))
            q = q / n;
        end

        eps = q(1:3);
        eta = q(4);
        I = sym(eye(3));

        C = (2*eta^2 - 1)*I + 2*(eps*eps.') - 2*eta*skew(eps);
    else
        q = q(:);
        n = norm(q);
        if n > 0
            q = q / n;
        end

        eps = q(1:3);
        eta = q(4);
        I = eye(3);

        C = (2*eta*eta - 1)*I + 2*(eps*eps.') - 2*eta*skew(eps);
    end
end