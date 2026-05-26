function q = DCM_to_quat(C)
% DCM_to_quat
% Converts DCM to quaternion q = [e1; e2; e3; eta]
% Matches quat_to_DCM convention:
% C = (2*eta^2 - 1)I + 2*eps*eps' - 2*eta*skew(eps)

C = reshape(double(C),3,3);

tr = trace(C);

if tr > 0
    S = 2*sqrt(max(tr + 1, 0));
    
    eta = 0.25*S;
    e1  = (C(2,3) - C(3,2))/S;
    e2  = (C(3,1) - C(1,3))/S;
    e3  = (C(1,2) - C(2,1))/S;

else
    if C(1,1) > C(2,2) && C(1,1) > C(3,3)
        S = 2*sqrt(max(1 + C(1,1) - C(2,2) - C(3,3), 0));
        
        e1  = 0.25*S;
        e2  = (C(1,2) + C(2,1))/S;
        e3  = (C(1,3) + C(3,1))/S;
        eta = (C(2,3) - C(3,2))/S;

    elseif C(2,2) > C(3,3)
        S = 2*sqrt(max(1 + C(2,2) - C(1,1) - C(3,3), 0));
        
        e1  = (C(1,2) + C(2,1))/S;
        e2  = 0.25*S;
        e3  = (C(2,3) + C(3,2))/S;
        eta = (C(3,1) - C(1,3))/S;

    else
        S = 2*sqrt(max(1 + C(3,3) - C(1,1) - C(2,2), 0));
        
        e1  = (C(1,3) + C(3,1))/S;
        e2  = (C(2,3) + C(3,2))/S;
        e3  = 0.25*S;
        eta = (C(1,2) - C(2,1))/S;
    end
end

q = [e1; e2; e3; eta];

% Normalize
n = norm(q);
if n > 0
    q = q/n;
else
    q = [0;0;0;1];
end



end