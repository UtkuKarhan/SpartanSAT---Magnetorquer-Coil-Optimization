function coil = coil_parameters()
coil.N = [200; 200; 200]; % turns
coil.A = [0.01; 0.01; 0.01]; % m^2 
coil.R = [4.4; 4.4; 4.4]; % ohms
coil.I_max = [0.15; 0.15; 0.15]; % amps
coil.V_max = [5; 5; 5]; % volts
% Dipole limits
coil.m_max = coil.N .* coil.A .* coil.I_max;
end