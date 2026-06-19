% Magnetic torque produced by magnetorquer
function tau = magnetorquer_torque(m, B)
tau = cross(m, B);
end