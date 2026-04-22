function C = Cx(a)
% Cx
% Rotation matrix about x-axis
%
% Inputs:
%   a        : rotation angle
%   symbolic : true for symbolic output, false for numeric

    c = cos(a);
    s = sin(a);
    C = [1, 0, 0;
         0, c, s;
         0,-s, c];
end