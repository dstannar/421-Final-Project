function r2 = reframe_vector(r1, F1, F2, DCM)
% reframe_vector
% Returns vector r1 (given in F1) as seen from F2 (r2)
%
% Inputs:
%   r1  : 3x1 vector
%   F1  : 3x3 first frame vectrix
%   F2  : 3x3 second frame vectrix
%   DCM : 3x3 direction cosine matrix
%
% You must enter either:
%   - DCM
% or
%   - F1 and F2

    if nargin < 4
        DCM = [];
    end
    if nargin < 3
        F2 = [];
    end
    if nargin < 2
        F1 = [];
    end

    if ~isempty(DCM)
        r2 = DCM * r1;
    elseif ~isempty(F1) && ~isempty(F2)
        C21 = F2.' * F1; % orthonormal rotation matrix
        r2 = C21 * r1;
    else
        error('You must enter either F1 and F2 or the DCM');
    end
end