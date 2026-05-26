function q3 = quatMult(q2, q1)
% Scalar-last convention:
%   q = [epsilon; eta]
% Inputs:
%   q2 = [epsilon2; eta2]
%   q1 = [epsilon1; eta1]
%
% Output:
%   q3 = [epsilon3; eta3]
%
% Book formula:
%   epsilon3 = eta2*epsilon1 + eta1*epsilon2 + epsilon1 x epsilon2
%   eta3     = eta1*eta2 - epsilon2' * epsilon1

q2 = reshape(q2,4,1);
q1 = reshape(q1,4,1);

eps2 = q2(1:3);
eta2 = q2(4);

eps1 = q1(1:3);
eta1 = q1(4);

% ross product: eps1 x eps2
eps1_cross_eps2 = [eps1(2)*eps2(3) - eps1(3)*eps2(2);
                   eps1(3)*eps2(1) - eps1(1)*eps2(3);
                   eps1(1)*eps2(2) - eps1(2)*eps2(1)];

eps3 = eta2*eps1 + eta1*eps2 + eps1_cross_eps2;
eta3 = eta1*eta2 - eps2.'*eps1;

q3 = [eps3; eta3];

end