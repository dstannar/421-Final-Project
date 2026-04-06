function a_x = skew(a)
% skew
% Return the skew-symmetric (cross-product) matrix a_x such that:
%   a_x * b = a x b

    a = a(:);

    ax = a(1);
    ay = a(2);
    az = a(3);

    a_x = [  0, -az,  ay;
            az,   0, -ax;
           -ay,  ax,   0];
end