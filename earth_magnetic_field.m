% IGRF-based Earth magnetic field in ECI frame [T]
function B_eci = earth_magnetic_field(r_eci, t)
epochUTC = [2026 1 1 0 0 0];   % [Y M D h m s]
% Current UTC as numeric date vector
utc_dt = datetime(epochUTC) + seconds(t);
utc_vec = datevec(utc_dt); 
% Convert ECI -> LLA
lla = eci2lla(r_eci.', utc_vec);   % [lat lon alt] = [deg deg m]
lat = lla(1);
lon = lla(2);
h   = lla(3);
dyear = decimal_year_local(utc_dt);
% XYZ = [North East Down]
B_ned_nT = igrfmagm(h, lat, lon, dyear, 14);
B_ned = B_ned_nT(:) * 1e-9;    % nT -> T
C_ned2ecef = ned2ecef_dcm_local(lat, lon);
B_ecef = C_ned2ecef * B_ned;
C_eci2ecef = dcmeci2ecef('IAU-2000/2006', utc_vec);
C_eci2ecef = C_eci2ecef(:,:,1);
B_eci = C_eci2ecef.' * B_ecef;
end
function dyear = decimal_year_local(utc_dt)
y0 = datetime(year(utc_dt),1,1,0,0,0);
y1 = datetime(year(utc_dt)+1,1,1,0,0,0);
dyear = year(utc_dt) + seconds(utc_dt - y0) / seconds(y1 - y0);
end


function C_ned2ecef = ned2ecef_dcm_local(lat_deg, lon_deg)

lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);
sLat = sin(lat);
cLat = cos(lat);
sLon = sin(lon);
cLon = cos(lon);

C_ned2ecef = [ -sLat*cLon, -sLon, -cLat*cLon;
                -sLat*sLon,  cLon, -cLat*sLon;
                 cLat,       0,    -sLat      ];
end