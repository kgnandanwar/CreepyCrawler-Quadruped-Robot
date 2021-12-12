clc
close all
clear all

syms P P1 P2
syms vx vy v_sq
syms K K1 K2
syms L
syms m1 m2 I1 I2 d1 d2 l1 l2 q th1 th2 dq dth1 dth2 ddq ddth1 ddth2
syms g real
syms dL_dth1 dL_dth2 DL_ddth1 DL_ddth2 dDL_dtddth2 dDL_dtddth1
syms u1 u2

P1 = m1*g*d1*cos(th1);
P2 = m2*g*(l1*cos(th1) + d2*cos(th1+th2));
P = P1 + P2;

vx = -l1*sin(th1)*dth1 - d2*sin(th1+th2)*(dth1+dth2);
vy = l1*cos(th1)*dth1 + d2*cos(th1+th2)*(dth1+dth2);

v_sq = simplify(vx^2 + vy^2);

K1 = 0.5*m1*d1*d1*(dth1^2) + 0.5*I1*(dth1^2);
K2 = 0.5*m2*v_sq + 0.5*I2*(dth1+dth2)^2;
K = simplify(K1 + K2);

L = simplify(K - P);


dL_dth1 = jacobian(L,(th1));
dL_dth2 = jacobian(L,(th2));
DL_ddth1 = jacobian(L,(dth1));
DL_ddth2 = jacobian(L,(dth2));

dDL_dtddth1 = jacobian(DL_ddth1,[th1; dth1; th2; dth2]) * [dth1; ddth1; dth2; ddth2];
dDL_dtddth2 = jacobian(DL_ddth2,[th1; dth1; th2; dth2]) * [dth1; ddth1; dth2; ddth2];

u1 = simplify(dDL_dtddth1 - dL_dth1)
u2 = simplify(dDL_dtddth2 - dL_dth2)

X = sym('X', [4,1]);
X(1) = th1;
X(2) = th2;
X(3) = dth1;
X(4) = dth2;

solution = solve([u1==0 u2 ==0],[ddth1 ddth2]);
disp(simplify(solution.ddth1))
disp(simplify(solution.ddth2))

dX = sym('dX', [4,1]);
dX(1) = X(3);
dX(2) = X(4);
dX(3) = simplify(solution.ddth1);
dX(4) = simplify(solution.ddth2);
% disp(dX(3))
% disp(dX(4))

th1 = deg2rad(30);
th2 = deg2rad(45);
[t, y] = ode45(@myode_swing, [0,10],[th1,th2,0,0]);
plot(t,y(:,1));
figure
plot(t,y(:,2));
figure
plot(t,y(:,3));
figure
plot(t,y(:,4));

function dX = myode_swing(t,X)
    m1 = 1;
    m2 = 1; 
    l1 = 1; 
    l2 = 1;
    d1 = 0.45;
    d2 = 0.45;
    I1 = 0.084;
    I2 = 0.084;
    g = 9.81;
    
    dX = zeros(4,1);
    X = num2cell(X);
    [th1, th2, dth1, dth2] = deal(X{:});
    if abs(th1)>2*pi
        th1 = mod(th1,2*pi);
    end
    
    if abs(th2)>2*pi
        th2 = mod(th2,2*pi);
    end
       
    dX(1) = dth1;
    dX(2) = dth2;

    
    
    dX(3) = (d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
    dX(4) = -(d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
 
%     dX(3) = (g*m2*d2*sin(th1 + th2) + g*l1*m2*sin(th1) + g*m1*d1*sin(th1))/(m2*l1^2 + 2*m2*cos(th2)*l1*d2 + m1*d1^2 + m2*d2^2 + I1 + I2);
%     dX(4) = (- l1*m2*d2*sin(th2)*dth1^2 + g*m2*d2*sin(th1 + th2))/(m2*d2^2 + I2);

end  
    
function dX = myode(t,X)
    m1 = 1;
    m2 = 1; 
    l1 = 1; 
    l2 = 1;
    d1 = 0.45;
    d2 = 0.45;
    I1 = 0.084;
    I2 = 0.084;
    g = 9.81;
    
    dX = zeros(4,1);
    X = num2cell(X);
    [th1, th2, dth1, dth2] = deal(X{:});
    if abs(th1)>2*pi
        th1 = mod(th1,2*pi);
    end
    
    if abs(th2)>2*pi
        th2 = mod(th2,2*pi);
    end
       
    dX(1) = dth1;
    dX(2) = dth2;
    dX(3) = (d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
    dX(4) = -(d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
 
%     dX(3) = (g*m2*d2*sin(th1 + th2) + g*l1*m2*sin(th1) + g*m1*d1*sin(th1))/(m2*l1^2 + 2*m2*cos(th2)*l1*d2 + m1*d1^2 + m2*d2^2 + I1 + I2);
%     dX(4) = (- l1*m2*d2*sin(th2)*dth1^2 + g*m2*d2*sin(th1 + th2))/(m2*d2^2 + I2);

end  

