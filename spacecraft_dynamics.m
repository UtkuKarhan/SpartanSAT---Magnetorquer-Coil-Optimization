function xdot = spacecraft_dynamics(t, x, J, k, coil)
q = x(1:4);
omega = x(5:7);
q = q / norm(q);
% Orbit position in ECI
r_eci = orbit_position(t);
% Earth magnetic field in Earth centered inertial frame
B_eci = earth_magnetic_field(r_eci, t);
% Rotate into body frame
R = quat_to_dcm(q);
B_body = R' * B_eci;
m_cmd = bdot_controller(omega, B_body, k, coil.m_max);
% m_cmd = [0.1; 0; 0];   % fixed x-axis dipole [A*m^2] (open-loop)
% Coil actuator producing the actual dipole
[m_actual, ~, ~, ~] = coil_actuator_model(m_cmd, coil);
% Torque
tau = cross(m_actual, B_body);
% Angular dynamics
omega_dot = J \ (tau - cross(omega, J*omega));
% Quaternions
Omega = [0        -omega(1) -omega(2) -omega(3);
         omega(1)  0         omega(3) -omega(2);
         omega(2) -omega(3)  0         omega(1);
         omega(3)  omega(2) -omega(1)  0       ];
q_dot = 0.5 * Omega * q;
xdot = [q_dot; omega_dot];
end