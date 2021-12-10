function T = twist2ht(S,theta)
    % your code here
    v = S(4:6);
    omega = S(1:3);
    T = [axisangle2rot(omega, theta) (eye(3)*theta + ((1-cos(theta)) * skew(omega)) + (theta - sin(theta))*skew(omega)^2)*v; 0 0 0 1];
    % If needed, you can calculate a rotation matrix with:
    % R = axisangle2rot(omega,theta);
end
function s = skew(a)
    % Writing a function to create a skew symmetric matrix
    s = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end
