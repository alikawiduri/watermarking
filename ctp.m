function [r, theta_rad, theta_deg] = ctp(x1, x2)
    % Menghitung radius (r)
    r = sqrt(x1^2 + x2^2);
    
    % Menghitung sudut (theta) dalam radian
    theta_rad = atan2(x2, x1);
    
    % Konversi sudut ke derajat
    theta_deg = rad2deg(theta_rad);
end