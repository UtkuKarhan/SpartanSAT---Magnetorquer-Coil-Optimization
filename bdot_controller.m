function m_cmd = bdot_controller(omega, B, k, ~)
m_cmd = k * cross(omega, B);
end