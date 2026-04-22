function C = Cz(a, symbolic)
% Cz
% Rotation matrix about z-axis
%
% Inputs:
%   a        : rotation angle
%   symbolic : true for symbolic output, false for numeric

    if symbolic
        c = cos(a);
        s = sin(a);
        C = [ c, s, 0;
             -s, c, 0;
              0, 0, 1];
    else
        c = cos(a);
        s = sin(a);
        C = [ c, s, 0;
             -s, c, 0;
              0, 0, 1];
    end
end