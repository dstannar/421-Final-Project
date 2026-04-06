function quat = DCM_to_quat(C)
% DCM_to_quat
% DCM to quaternion
% Returns q = [e1; e2; e3; eta], normalized.
%
% Supports numeric or symbolic input.

    tol = 1e-8;

    if isa(C, 'sym')
        tr = trace(C);
        eta = sqrt(tr + 1)/2;

        if abs(double(vpa(eta))) > tol
            e1 = (C(3,2) - C(2,3)) / (4*eta);
            e2 = (C(1,3) - C(3,1)) / (4*eta);
            e3 = (C(2,1) - C(1,2)) / (4*eta);
        else
            e1 = sqrt((C(1,1) + 1)/2);
            e2 = sqrt((C(2,2) + 1)/2);
            e3 = sqrt((C(3,3) + 1)/2);

            % arbitrarily choose e1 positive
            e1 = abs(e1);
            e2 = sign(C(1,2)) * abs(e2);
            e3 = sign(C(1,3)) * abs(e3);
        end

        quat = [e1; e2; e3; eta];
        quat = quat / norm(quat);

    else
        C = reshape(double(C), 3, 3);
        tr = trace(C);
        eta = 0.5 * sqrt(tr + 1.0);

        if abs(eta) > tol
            e1 = (C(2,3) - C(3,2)) / (4.0*eta);
            e2 = (C(3,1) - C(1,3)) / (4.0*eta);
            e3 = (C(2,1) - C(1,2)) / (4.0*eta);
        else
            e1 = sqrt((C(1,1) + 1.0)/2.0);
            e2 = sqrt((C(2,2) + 1.0)/2.0);
            e3 = sqrt((C(3,3) + 1.0)/2.0);

            % arbitrarily choose e1 positive
            e1 = abs(e1);
            e2 = sign(C(1,2)) * abs(e2);
            e3 = sign(C(1,3)) * abs(e3);
        end

        quat = [e1; e2; e3; eta];
        quat = quat / norm(quat);
    end
end