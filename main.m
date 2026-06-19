% PARAMETERS
m = 2.66;
% Inertia matrix [kg*m^2]
Jx = (1/12)*m*(0.1^2 + 0.2^2);
Jy = (1/12)*m*(0.1^2 + 0.2^2);
Jz = (1/12)*m*(0.1^2 + 0.1^2);
J = diag([Jx Jy Jz]);
% Initial attitude quaternion [q0 q1 q2 q3]^T
q0 = [1; 0; 0; 0]; % spacecraft starts aligned with the inertial frame (no rotation)
% Initial angular velocity [rad/s]
omega0 = deg2rad([8; -6; 5]); % The initial angular velocity was assumed to represent a post-deployment tumbling 
% Initial state vector
x0 = [q0; omega0];
% Simulation time span [s]
tspan = [0 1500];
% B-dot controller gain
k = 5e4;
% k = 0; 
coil = coil_parameters();
% RUN SIMULATION
[t, x] = ode45(@(t,x) spacecraft_dynamics(t, x, J, k, coil), tspan, x0);
% EXTRACT STATES
q_hist = x(:,1:4);
omega = x(:,5:7);
omega_mag = vecnorm(omega, 2, 2);
% STATE PLOTS
% Angular velocity components
figure;
plot(t, rad2deg(omega(:,1)), 'LineWidth', 1.2); hold on;
plot(t, rad2deg(omega(:,2)), 'LineWidth', 1.2);
plot(t, rad2deg(omega(:,3)), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Angular Velocity (deg/s)');
legend('\omega_x', '\omega_y', '\omega_z');
title('Angular Velocity Components');
grid on;
% Angular velocity magnitude
figure;
plot(t, rad2deg(omega_mag), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('|\omega| (deg/s)');
title('Angular Velocity Magnitude');
grid on;
% Rotational kinetic energy
E = zeros(length(t),1);
for i = 1:length(t)
    w = omega(i,:)';
    E(i) = 0.5 * w' * J * w;
end
figure;
plot(t, E, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Rotational Kinetic Energy (J)');
title('Rotational Kinetic Energy');
grid on;
% COIL / ACTUATOR LOGGING
Npts = length(t);
m_cmd_hist    = zeros(Npts, 3);
m_actual_hist = zeros(Npts, 3);
I_hist        = zeros(Npts, 3);
V_hist        = zeros(Npts, 3);
P_hist        = zeros(Npts, 3);
B_body_hist   = zeros(Npts, 3);
for i = 1:Npts
    q = q_hist(i,:)';
    w = omega(i,:)';
    % Orbit position in ECI
    r_eci = orbit_position(t(i));
    % Magnetic field in ECI
    B_eci = earth_magnetic_field(r_eci, t(i));
    % Convert magnetic field to body frame
    R = quat_to_dcm(q);
    B_body = R' * B_eci;
    B_body_hist(i,:) = B_body';
    % Controller-commanded dipole
    m_cmd = bdot_controller(w, B_body, k, coil.m_max);
    % m_cmd = [0.1; 0; 0];   % open-loop fixed x-axis dipole
    % Actual coil response
    [m_actual, I_cmd, V_cmd, P_cmd] = coil_actuator_model(m_cmd, coil);
    m_cmd_hist(i,:)    = m_cmd';
    m_actual_hist(i,:) = m_actual';
    I_hist(i,:)        = I_cmd';
    V_hist(i,:)        = V_cmd';
    P_hist(i,:)        = P_cmd';
end
P_total = sum(P_hist, 2);
% COIL / ACTUATOR PLOTS
% Commanded vs actual dipole moment
figure;
% X-axis dipole
plot(t, m_cmd_hist(:,1), 'k--', 'LineWidth', 2); hold on;
plot(t, m_actual_hist(:,1), 'r-', 'LineWidth', 1.5);
% Y-axis dipole
plot(t, m_cmd_hist(:,2), 'b--', 'LineWidth', 2);
plot(t, m_actual_hist(:,2), 'm-', 'LineWidth', 1.5);
% Z-axis dipole
plot(t, m_cmd_hist(:,3), 'g--', 'LineWidth', 2);
plot(t, m_actual_hist(:,3), 'c-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Dipole Moment (A·m^2)');
legend('m_{cmd,x}', 'm_{act,x}', ...
       'm_{cmd,y}', 'm_{act,y}', ...
       'm_{cmd,z}', 'm_{act,z}');
title('Commanded vs Actual Dipole Moment');
grid on;
% Coil current
figure;
plot(t, I_hist(:,1), 'LineWidth', 1.2); hold on;
plot(t, I_hist(:,2), 'LineWidth', 1.2);
plot(t, I_hist(:,3), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Current (A)');
legend('I_x', 'I_y', 'I_z');
title('Coil Current');
grid on;
% Coil voltage
figure;
plot(t, V_hist(:,1), 'LineWidth', 1.2); hold on;
plot(t, V_hist(:,2), 'LineWidth', 1.2);
plot(t, V_hist(:,3), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_x', 'V_y', 'V_z');
title('Coil Voltage');
grid on;
% Coil power per axis
figure;
plot(t, P_hist(:,1), 'LineWidth', 1.2); hold on;
plot(t, P_hist(:,2), 'LineWidth', 1.2);
plot(t, P_hist(:,3), 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Power (W)');
legend('P_x', 'P_y', 'P_z');
title('Coil Power by Axis');
grid on;
% Total coil power
figure;
plot(t, P_total, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Total Coil Power (W)');
title('Total Coil Power Consumption');
grid on;
% PERFORMANCE METRICS
omega0_deg = rad2deg(omega_mag(1));
omegaf_deg = rad2deg(omega_mag(end));
threshold_deg = 0.5;
idx_detumble = find(rad2deg(omega_mag) < threshold_deg, 1, 'first');
if isempty(idx_detumble)
    detumble_time = NaN;
else
    detumble_time = t(idx_detumble);
end
peak_power = max(P_total);
peak_current = max(abs(I_hist), [], 1);
peak_voltage = max(abs(V_hist), [], 1);
peak_dipole_actual = max(abs(m_actual_hist), [], 1);
fprintf('Initial angular speed: %.3f deg/s\n', omega0_deg);
fprintf('Final angular speed:   %.3f deg/s\n', omegaf_deg);
if isnan(detumble_time)
    fprintf('Detumble threshold of %.2f deg/s was not reached.\n', threshold_deg);
else
    fprintf('Detumble time to %.2f deg/s: %.3f s\n', threshold_deg, detumble_time);
end
fprintf('Peak total coil power: %.4f W\n', peak_power);
fprintf('Peak current [Ix Iy Iz]: [%.4f %.4f %.4f] A\n', ...
    peak_current(1), peak_current(2), peak_current(3));
fprintf('Peak voltage [Vx Vy Vz]: [%.4f %.4f %.4f] V\n', ...
    peak_voltage(1), peak_voltage(2), peak_voltage(3));
fprintf('Peak actual dipole [mx my mz]: [%.4f %.4f %.4f] A*m^2\n', ...
    peak_dipole_actual(1), peak_dipole_actual(2), peak_dipole_actual(3));
m_error = m_cmd_hist - m_actual_hist;
figure;
plot(t, m_error(:,1), 'LineWidth', 1.5); hold on;
plot(t, m_error(:,2), 'LineWidth', 1.5);
plot(t, m_error(:,3), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Dipole Error (A·m^2)');
legend('e_x','e_y','e_z');
title('Commanded - Actual Dipole Error');
grid on;
% ORBIT CHECK FOR THE ASSUMPTIONS VERIFICATION PART
r_hist = zeros(length(t),3);
for i = 1:length(t)
    r_hist(i,:) = orbit_position(t(i))';
end
figure;
plot3(r_hist(:,1), r_hist(:,2), r_hist(:,3), 'LineWidth', 1.5);
xlabel('x_{ECI} (m)');
ylabel('y_{ECI} (m)');
zlabel('z_{ECI} (m)');
title('Spacecraft Orbit in ECI');
grid on;
axis equal;