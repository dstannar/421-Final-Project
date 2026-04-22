function C = Cy(a, symbolic)
% Cy
% Rotation matrix about y-axis
%
% Inputs:
%   a        : rotation angle
%   symbolic : true for symbolic output, false for numeric

    if symbolic
        c = cos(a);
        s = sin(a);
        C = [ c, 0,-s;
              0, 1, 0;
              s, 0, c];
    else
        c = cos(a);
        s = sin(a);
        C = [ c, 0,-s;
              0, 1, 0;
              s, 0, c];
    end
end