function Tb_b = WMM(m_b, r_ECI, q_b_ECI, JD_0, t)


% DGRF 2020 (given)
% a = 6371.2; %km
% g_11 = -1450.9; %nT
% h_11 = 4652.5; %nT
% g_10 = -29404.8; %nT
% m_ECEF = a^3*[g_11; h_11; g_10];

% converting given values to m and T
a = 6371.2 * 10^3;
g_11 = -1450.9 * 10^-9; 
h_11 = 4652.5 * 10^-9; 
g_10 = -29404.8 * 10^-9; 
m_ECEF = a^3*[g_11; h_11; g_10];

% local time is noon at initial given state, assuming t is in seconds
JD = JD_0 + t / 86400;

% first convert r from eci to ecef 

C_ecef_eci = eci2ecef(JD, t);

r_ECEF = C_ecef_eci * r_ECI;

r = norm(r_ecef);

B_ecef = (3 * (m_ECEF.' * r_ECEF) * r_ECEF - r^2 * m_ECEF) / r^5;

C_b_ECI = quat2C(q_b_ECI);

C_b_ecef = C_b_ECI * C_ecef_eci'; % doing successive rotions and transposing to get C_b_eci * C_eci_ecef

B_b = C_b_ecef * B_ecef;

Tb_b = vcross(m_b) * B_b;



    function C_ecef_eci = eci2ecef(JD) 
        T0 = (JD - 2451545) / 36525; 
        theta_G = wrapTo360(100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - 2.58*10^-8 * T0^3); % in degs 
    
        C_ecef_eci = [cos(theta_G) sin(theta_G) 0
            -sin(theta_G) cos(theta_G) 0
            0 0 1]; % negative rotation about z axis 

    end 

        function vCross = vcross(v)
            vCross = [ 0,    -v(3),  v(2)
                      v(3),    0,   -v(1)
                     -v(2),   v(1),   0];
        end
        

        function C21 = quat2C(quaternion)
           
            epsilon = quaternion(1:3);
            eta = quaternion(4);
            I = [1 0 0
                0 1 0
                0 0 1];
            C21 = (2 * eta^2 - 1) * I + 2 * (epsilon * epsilon') - 2 * eta * vcross(epsilon);

        
        end

end

