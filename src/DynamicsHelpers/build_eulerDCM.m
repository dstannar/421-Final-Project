function DCM = build_eulerDCM(first_rot, second_rot, third_rot, a1, a2, a3, degs)
% build_eulerDCM
% Builds direction cosine matrix for any 3-axis rotation sequence
% Based off fundamentals of d&c, ch. 1.3.3
%
% Note:
%   1 = x-axis
%   2 = y-axis
%   3 = z-axis
%
% You must enter the rotation angles in the order they occur.
%
% Inputs:
%   first_rot  : first rotation axis number
%   second_rot : second rotation axis number
%   third_rot  : third rotation axis number
%   a1         : first rotation angle
%   a2         : second rotation angle
%   a3         : third rotation angle
%   degs       : true if input angles are in degrees
%
% Output:
%   DCM        : 3x3 direction cosine matrix

    if nargin < 7
        degs = false;
    end

    if nargin < 6
        a3 = [];
    end
    if nargin < 5
        a2 = [];
    end
    if nargin < 4
        a1 = [];
    end

    if isempty(a1) && isempty(a2) && isempty(a3)
        symbolic = true;

        angle_names = {'phi', 'theta', 'psi'};
        a1 = sym(angle_names{first_rot}, 'real');
        a2 = sym(angle_names{second_rot}, 'real');
        a3 = sym(angle_names{third_rot}, 'real');
    else
        symbolic = false;
    end

    if degs && ~symbolic
        a1 = deg2rad(a1);
        a2 = deg2rad(a2);
        a3 = deg2rad(a3);
    end

    switch first_rot
        case 1
            matrix1 = Cx(a1);
        case 2
            matrix1 = Cy(a1);
        case 3
            matrix1 = Cz(a1);
        otherwise
            error('first_rot must be 1, 2, or 3');
    end

    switch second_rot
        case 1
            matrix2 = Cx(a2);
        case 2
            matrix2 = Cy(a2);
        case 3
            matrix2 = Cz(a2);
        otherwise
            error('second_rot must be 1, 2, or 3');
    end

    switch third_rot
        case 1
            matrix3 = Cx(a3);
        case 2
            matrix3 = Cy(a3);
        case 3
            matrix3 = Cz(a3);
        otherwise
            error('third_rot must be 1, 2, or 3');
    end

    DCM = matrix3 * matrix2 * matrix1;
end