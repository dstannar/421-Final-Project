function q_conj = quatConjugate(q)
% quatConjugate Quaternion conjugate using scalar-last convention.
%
% Input:
%   q = [epsilon; eta] = [q1; q2; q3; q4]
%
% Output:
%   q_conj = [-epsilon; eta]

q = reshape(q,4,1);

q_conj = [-q(1:3); q(4)];

end