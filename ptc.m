function [x1, x2] = ptc(r, theta_deg)
    % Konversi sudut dari derajat ke radian
    theta_rad = deg2rad(theta_deg);
    
    % Menghitung koordinat kartesian
    x1 = r * cos(theta_rad);
    x2 = r * sin(theta_rad);
end
