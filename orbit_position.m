function r_eci = orbit_position(t)
% Earth constants
mu = 3.986004418e14;   % Earth's gravitational parameter [m^3/s^2]
Re = 6371e3;           % Earth radius [m]
% Mission orbital parameters
h_perigee = 500e3;     % perigee altitude [m]
h_apogee  = 500e3;     % apogee altitude [m]
i         = deg2rad(51.6);   % inclination [rad]
RAAN  = deg2rad(0);    % right ascension of ascending node [rad] 
omega = deg2rad(0);    % argument of perigee [rad] 
nu0   = deg2rad(0);    % initial true anomaly [rad] 
% Convert altitudes to orbital radii
rp = Re + h_perigee;   % perigee radius [m]
ra = Re + h_apogee;    % apogee radius [m]
% Orbital elements
a = (rp + ra)/2; % semi-major axis [m]
e = (ra - rp)/(ra + rp); % eccentricity
p = a * (1 - e^2); % semi-latus rectum [m]
% Mean motion
n = sqrt(mu / a^3); % mean motion [rad/s]
% Propagate mean anomaly
M = n * t;
% Kepler's equation: M = E - e*sin(E)
E = kepler_solver(M, e);
% True anomaly from eccentric anomaly
nu = 2 * atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
nu = nu0 + nu;
% Perifocal position
r_pf = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
% Rotation from perifocal frame to ECI
R3_RAAN = [ cos(RAAN) -sin(RAAN) 0;
            sin(RAAN)  cos(RAAN) 0;
            0          0         1];
R1_i = [ 1  0       0;
         0  cos(i) -sin(i);
         0  sin(i)  cos(i)];
R3_omega = [ cos(omega) -sin(omega) 0;
             sin(omega)  cos(omega) 0;
             0           0          1];
Q_pX = R3_RAAN * R1_i * R3_omega;
% ECI position
r_eci = Q_pX * r_pf;
end
% KEPLER SOLVER
function E = kepler_solver(M, e)
% Newton-Raphson
E = M; 
for k = 1:20
    f = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    E = E - f/fp;
end
end