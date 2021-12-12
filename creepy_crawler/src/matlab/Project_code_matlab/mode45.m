
function dy = mode45(t,y)

m1 = 1;
m2 = 1;
l1 = 1; 
l2 = 1;
r1 = 0.45;
r2 = 0.45;
I1 = 0.084;
I2 = 0.084;
g = 9.81;


y = num2cell(y);
[q1,q2,q1dot,q2dot]=deal(y{:});

if abs(q1) > 2*pi
    q1 = mod(q1,2*pi);
end

if abs(q2) > 2*pi
    q2 = mod(q2,2*pi);
end

dy = zeros(4,1);
dy(1) = q1dot;
dy(2) = q2dot;
dy(3) = (I2*g*l1*m2*sin(q2) - m2^2*q1dot*q2dot*r2^4*sin(2*q2) - 0.5000*I2*m2*q1dot^2*r2^2*sin(2*q2) - I2*m2*q1dot*q2dot*r2^2*sin(2*q2) + 2*l1*m2^2*q1dot*q2dot*r2^3*cos(q2) - I2*l1*m2*q1dot^2*r2*cos(q2) - 4*l1*m2^2*q1dot*q2dot*r2^3*cos(q1)^2*cos(q2) - 2*l1^2*m2^2*q1dot^2*r2^2*cos(q1)*cos(q2)^2*sin(q1) + 2*I2*l1*m2*q1dot*q2dot*r2*cos(q2) + 2*I2*l1*m2*q1dot^2*r2*cos(q1)^2*cos(q2) + 4*l1*m2^2*q1dot^2*r2^3*cos(q1)*sin(q1)*sin(q2) + 2*l1*m2^2*q2dot^2*r2^3*cos(q1)*sin(q1)*sin(q2) + 4*l1^2*m2^2*q1dot^2*r2^2*cos(q1)^3*cos(q2)^2*sin(q1) - 2*l1*m2^2*q1dot^2*r2^3*cos(q1)*cos(q2)^2*sin(q1)*sin(q2) + 2*g*l1^2*m2^2*r2*cos(q1)*cos(q2)*sin(q1)*sin(q2) + 4*I2*l1*m2*q1dot^2*r2*cos(q1)*sin(q1)*sin(q2) + 2*I2*l1*m2*q2dot^2*r2*cos(q1)*sin(q1)*sin(q2) - 4*I2*l1*m2*q1dot*q2dot*r2*cos(q1)^2*cos(q2))/(I1*I2 + m2^2*r2^4 - m2^2*r2^4*cos(q2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + 2*I2*m2*r2^2 - I2*m2*r2^2*cos(q2)^2 - 2*l1*m2^2*r2^3*sin(q2) + m1*m2*r1^2*r2^2 - 4*l1^2*m2^2*r2^2*cos(q1)^2*cos(q2)^2 + 4*l1^2*m2^2*r2^2*cos(q1)^4*cos(q2)^2 + 4*l1*m2^2*r2^3*cos(q1)^2*sin(q2) - 2*I2*l1*m2*r2*sin(q2) + 4*I2*l1*m2*r2*cos(q1)^2*sin(q2) - 4*I2*l1*m2*r2*cos(q1)*cos(q2)*sin(q1));
dy(4) = -(g*l1^3*m2^2*sin(q2) - 0.5000*m2^2*q1dot^2*r2^4*sin(2*q2) - 2*g*l1^2*m2^2*r2 + 4*g*l1^2*m2^2*r2*cos(q1)^2 + 2*g*l1^2*m2^2*r2*cos(q2)^2 + l1*m2^2*q1dot^2*r2^3*cos(q2) - l1^3*m2^2*q1dot^2*r2*cos(q2) + g*l1*m2^2*r2^2*sin(q2)^3 - l1*m2^2*q1dot^2*r2^3*cos(q2)^3 + I1*g*l1*m2*sin(q2) + I2*g*l1*m2*sin(q2) + 0.5000*l1^2*m2^2*q1dot^2*r2^2*sin(2*q2) + m2^2*q1dot^2*r2^4*cos(q2)^3*sin(q2) - 0.5000*I1*m2*q1dot^2*r2^2*sin(2*q2) - 0.5000*I2*m2*q1dot^2*r2^2*sin(2*q2) - I2*m2*q1dot*q2dot*r2^2*sin(2*q2) - 4*g*l1^2*m2^2*r2*cos(q1)^2*cos(q2)^2 - 2*l1*m2^2*q1dot^2*r2^3*cos(q1)^2*cos(q2) + 2*l1^3*m2^2*q1dot^2*r2*cos(q1)^2*cos(q2) - I1*l1*m2*q1dot^2*r2*cos(q2) - I2*l1*m2*q1dot^2*r2*cos(q2) + 2*l1*m2^2*q1dot^2*r2^3*cos(q1)^2*cos(q2)^3 + g*l1*m1*m2*r1^2*sin(q2) - 0.5000*m1*m2*q1dot^2*r1^2*r2^2*sin(2*q2) + 4*l1^2*m2^2*q2dot^2*r2^2*cos(q1)^2*cos(q2)*sin(q2) - 4*l1^2*m2^2*q2dot^2*r2^2*cos(q1)^4*cos(q2)*sin(q2) - l1*m1*m2*q1dot^2*r1^2*r2*cos(q2) + 2*I2*l1*m2*q1dot*q2dot*r2*cos(q2) + 2*I1*l1*m2*q1dot^2*r2*cos(q1)^2*cos(q2) + 2*I2*l1*m2*q1dot^2*r2*cos(q1)^2*cos(q2) + 4*I2*l1*m2*q1dot^2*r2*cos(q1)*sin(q1)*sin(q2) + 2*I2*l1*m2*q2dot^2*r2*cos(q1)*sin(q1)*sin(q2) + 4*l1^2*m2^2*q1dot*q2dot*r2^2*cos(q1)*cos(q2)^2*sin(q1) + 2*l1*m1*m2*q1dot^2*r1^2*r2*cos(q1)^2*cos(q2) - 4*I2*l1*m2*q1dot*q2dot*r2*cos(q1)^2*cos(q2) - 8*l1^2*m2^2*q1dot*q2dot*r2^2*cos(q1)^3*cos(q2)^2*sin(q1) - 4*l1*m2^2*q1dot*q2dot*r2^3*cos(q1)*cos(q2)^2*sin(q1)*sin(q2))/(I1*I2 + m2^2*r2^4 - m2^2*r2^4*cos(q2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 + 2*I2*m2*r2^2 - I2*m2*r2^2*cos(q2)^2 - 2*l1*m2^2*r2^3*sin(q2) + m1*m2*r1^2*r2^2 - 4*l1^2*m2^2*r2^2*cos(q1)^2*cos(q2)^2 + 4*l1^2*m2^2*r2^2*cos(q1)^4*cos(q2)^2 + 4*l1*m2^2*r2^3*cos(q1)^2*sin(q2) - 2*I2*l1*m2*r2*sin(q2) + 4*I2*l1*m2*r2*cos(q1)^2*sin(q2) - 4*I2*l1*m2*r2*cos(q1)*cos(q2)*sin(q1));

end

