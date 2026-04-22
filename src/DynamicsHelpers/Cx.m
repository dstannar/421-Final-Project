function C = Cx(a, symbolic)
% Cx
% Rotation matrix about x-axis
%
% Inputs:
%   a        : rotation angle
%   symbolic : true for symbolic output, false for numeric

    if symbolic
        c = cos(a);
        s = sin(a);
        C = [1, 0, 0;
             0, c, s;
             0,-s, c];
    else
        c = cos(a);
        s = sin(a);
        C = [1, 0, 0;
             0, c, s;
             0,-s, c];
    end
end