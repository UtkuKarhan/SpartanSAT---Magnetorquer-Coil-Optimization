function [m_actual, I_cmd, V_cmd, P_cmd] = coil_actuator_model(m_cmd, coil)
% Convert desired dipole to required current
I_req = m_cmd ./ (coil.N .* coil.A);
% Saturate by current limits
I_cmd = max(min(I_req, coil.I_max), -coil.I_max);
% Compute voltage needed
V_req = I_cmd .* coil.R;
% voltage limits
overV = abs(V_req) > coil.V_max;
I_cmd(overV) = sign(I_cmd(overV)) .* (coil.V_max(overV) ./ coil.R(overV));
% Actual voltage after saturation
V_cmd = I_cmd .* coil.R;
% Actual dipole produced
m_actual = coil.N .* coil.A .* I_cmd;
% Coil power
P_cmd = I_cmd.^2 .* coil.R;
end