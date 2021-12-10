function R = axisangle2rot(omega,theta)
    % Applying Rodrigues' Formula 
    R = eye(3) + (sin(theta) * skew(omega)) + ((1 - cos(theta)) * (skew(omega))^2);
end