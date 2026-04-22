function [C_LVLH_ECI] = LVLH_from_ECI(r_eci, v_eci)
% getLVLHFrame
% Build LVLH basis vectors and the DCM relating ECI to LVLH

    z_lvlh_eci = -r_eci / norm(r_eci);
    y_lvlh_eci = -cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
    x_lvlh_eci = cross(y_lvlh_eci, z_lvlh_eci);

    C_LVLH_ECI = [x_lvlh_eci.';
                  y_lvlh_eci.';
                  z_lvlh_eci.'];
end