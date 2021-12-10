clc
close all
clear all
sympref('FloatingPointOutput', true);

syms P P1 P2
syms vx vy v_sq
syms K K1 K2
syms L
% Consider d1 = r1 and d2 = r2
syms m1 m2 I1 I2 d1 d2 l1 l2 q th1 th2 dq dth1 dth2 ddq ddth1 ddth2
syms g real
syms dL_dth1 dL_dth2 DL_ddth1 DL_ddth2 dDL_dtddth2 dDL_dtddth1
syms u lambda tau1 tau2

m1 = 1;
m2 = 1; 
l1 = 1; 
l2 = 1;
d1 = 0.45;
d2 = 0.45;
I1 = 0.084;
I2 = 0.084;
g = 9.81;

Z(1) = dth1;
Z(2) = dth2;
Z(3) = simplify((I2*tau1 - I2*tau2 + d2^2*m2*tau1 - d2^2*m2*tau2 + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) - d2*l1*m2*tau2*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2));
Z(4) = simplify(-(I2*tau1 - I1*tau2 - I2*tau2 - d1^2*m1*tau2 + d2^2*m2*tau1 - d2^2*m2*tau2 - l1^2*m2*tau2 + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + d2*l1*m2*tau1*cos(th2) - 2*d2*l1*m2*tau2*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2));

% Equilibrium Points: Tau1=Tau2=0
dX(1) = dth1;
dX(2) = dth2;
dX(3) = simplify(-(d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2));dX(4) = simplify(-(I2*tau1 - I1*tau2 - I2*tau2 - d1^2*m1*tau2 + d2^2*m2*tau1 - d2^2*m2*tau2 - l1^2*m2*tau2 + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + d2*l1*m2*tau1*cos(th2) - 2*d2*l1*m2*tau2*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2));
dX(4) = simplify(-(d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2));

sol = solve([dX(1) ==0, dX(2)==0, dX(3)==0, dX(4)==0], [dth1, dth2, th1, th2]);
dth1_Sol = sol.dth1;
dth2_Sol = sol.dth2;
th1_Sol = sol.th1;
th2_Sol = sol.th2;


Z = [Z(1); Z(2); Z(3); Z(4)];
X = [th1; th2; dth1; dth2];
u = [tau1; tau2];

% LINEARIZATION
A = jacobian(Z,X);
B = jacobian(Z, [u]);

A = subs(A,[th1, th2, dth1, dth2],[th1_Sol, th2_Sol, dth1_Sol, dth2_Sol]);
B = subs(B,[th1, th2, dth1, dth2],[th1_Sol, th2_Sol, dth1_Sol, dth2_Sol]);

A = double(A);
B = double(B);

% STABILITY
eigA = eig(A);

% CONTROLLABILITY
rankCO = rank(ctrb(A,B));

% STATE FEEDBACK DESIGN  FOR -1, -2, -1-i, -1+i
lambda = [-1,-2,-1-1i,-1+1i];
Kn = place(A,B,lambda);

% ODE
th1 = deg2rad(30);
th2 = deg2rad(45);
dth1 = 0;
dth2 = 0;
[t, y] = ode45(@myode, [0,10],[th1,th2,dth1,dth2]);
plot(t,y(:,1));
figure
plot(t,y(:,2));
figure
plot(t,y(:,3));
figure
plot(t,y(:,4));

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

    K1 = [23.5864  ,  5.8874 ,   5.1472,    2.6109];
    K2 = [5.8879  ,  4.9875 ,   1.5444  ,  0.9770];

    tau1 = -K1 * [th1; th2; dth1; dth2];
    tau2 = -K2 * [th1; th2; dth1; dth2];
       
    dX(1) = dth1;
    dX(2) = dth2;
    dX(3) = (I2*tau1 - I2*tau2 + d2^2*m2*tau1 - d2^2*m2*tau2 + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) - d2*l1*m2*tau2*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
    dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - d1^2*m1*tau2 + d2^2*m2*tau1 - d2^2*m2*tau2 - l1^2*m2*tau2 + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + d2*l1*m2*tau1*cos(th2) - 2*d2*l1*m2*tau2*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
end    
