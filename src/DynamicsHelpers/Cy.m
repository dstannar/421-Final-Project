function C = Cy(a)
% Cy
% Rotation matrix about y-axis
%
% Inputs:
%   a        : rotation angle
%   symbolic : true for symbolic output, false for numeric

    c = cos(a);
    s = sin(a);
    C = [ c, 0,-s;
          0, 1, 0;
          s, 0, c];
end