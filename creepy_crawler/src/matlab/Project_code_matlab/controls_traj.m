clc
close all
clear all

syms q k l m n t 
syms k1 l1 m1 n1 k2 l2 m2 n2
syms m1 m2 I1 I2 d1 d2 l1 l2 q th1 th2 dq dth1 dth2 ddq ddth1 ddth2
syms g real
syms K v
syms tau1 tau2

%------------------------------------------------------------------------------------------------------
% Kindly Comment other sections while running a particular code, say,
% comment part b,c,d,e,f while checking and running code for part a 
%------------------------------------------------------------------------------------------------------
% % Part a

% t0 = 0;
% tf = 10;
% th10 = deg2rad(180);
% th1f = 0;
% th20 = deg2rad(90);
% th2f = 0;
% dth10 = 0;
% dth1f = 0;
% dth20 = 0;
% dth2f = 0;
% 
% time = [t0 tf];
% th1 = [th10 th1f];
% th2 = [th20 th2f];
% th = [th1; th2];
% dth1 = [dth10 dth1f];
% dth2 = [dth20 dth2f];
% dth = [dth1; dth2];
% k = [k1; k2];
% l = [l1 ;l2];
% m = [m1 ;m2];
% n = [n1 ;n2];
% 
% 
% eq = k + l*t + m*t^2 + n*t^3;
% q = subs(eq, t, time);
% q_dot = subs(diff(eq),t,time);
% sol = solve([q==th q_dot==dth], [k,l,m,n]);
% ans_a = subs(eq, sol)                              % final answer (a)
%------------------------------------------------------------------------------------------------------
% Part b

% th = [th1; th2];
% dth = [dth1; dth2];
% ddth = [ddth1; ddth2];
m1 = 1;
m2 = 1; 
l1 = 1; 
l2 = 1;
d1 = 0.45;
d2 = 0.45;
I1 = 0.084;
I2 = 0.084;
g = 9.81;
u1 = ddth1*(m1*d1^2 + I1 + I2 + 0.5000*m2*(2*d2^2 + 4*cos(th2)*d2*l1 + 2*l1^2)) + ddth2*(I2 + 0.5000*m2*(2*d2^2 + 2*l1*cos(th2)*d2)) - g*(m2*(d2*sin(th1 + th2) + l1*sin(th1)) - d1*m1*sin(th1))  - 0.5000*dth2*m2*(4*d2*dth1*l1*sin(th2) + 2*d2*dth2*l1*sin(th2));
u2 = ddth1*(I2 + d2^2*m2 + d2*l1*m2*cos(th2)) + ddth2*(I2 + d2^2*m2) + th1*(d2*dth1*l1*m2*sin(th2)) - d2*g*m2*sin(th1 + th2) ;
% tau = [tau1;tau2];
% u = [u1;u2];
% % u = M*ddth + C*dth + G == tau;                       %Equation
G = [- g*(m2*(d2*sin(th1 + th2) + l1*sin(th1)) - d1*m1*sin(th1)); - d2*g*m2*sin(th1 + th2) ];
C = [0, 0.5000*m2*(4*d2*dth1*l1*sin(th2) + 2*d2*dth2*l1*sin(th2)) ; (d2*dth1*l1*m2*sin(th2)), 0];
M = [(m1*d1^2 + I1 + I2 + 0.5000*m2*(2*d2^2 + 4*cos(th2)*d2*l1 + 2*l1^2)), (I2 + 0.5000*m2*(2*d2^2 + 2*l1*cos(th2)*d2)); (I2 + d2^2*m2 + d2*l1*m2*cos(th2)), (I2 + d2^2*m2)];
% %-----------------------------------------------------------------------------------------------------------------------------------------------------
% % Part c
% % Tau = M*v + C*dth + G

%-------------------------------------------------------------------------------------------------------------------

% part d
% ODE

th1 = deg2rad(200);
th2 = deg2rad(125);
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

A=[ 0  0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
B=[0 0; 0 0 ; 1 0; 0 1];

lambda=[-10,-100,-10,-100];
K=place(A,B,lambda);

     q1 = 0.0063*t^3 - 0.0942*t^2 + 3.1416;
     q2 = 0.0031*t^3 - 0.0471*t^2 + 1.5708;

     dq1 = 3*0.0063*t^2 - 2*0.0942*t;
     dq2 = 3*0.0031*t^2 - 2*0.0471*t;

     ddq1  = 6*0.0063*t - 2*0.0942;
     ddq2 = 6*0.0031*t - 2*0.0471;

     G = [- g*(m2*(d2*sin(th1 + th2) + l1*sin(th1)) - d1*m1*sin(th1)); - d2*g*m2*sin(th1 + th2) ];
     C = [0, 0.5000*m2*(4*d2*dth1*l1*sin(th2)) + 0.5000*m2*(2*d2*dth2*l1*sin(th2)) ; (d2*dth1*l1*m2*sin(th2)), 0];
     M = [(m1*d1^2 + I1 + I2 + 0.5000*m2*(2*d2^2 + 4*cos(th2)*d2*l1 + 2*l1^2)), (I2 + 0.5000*m2*(2*d2^2 + 2*l1*cos(th2)*d2)); (I2 + d2^2*m2 + d2*l1*m2*cos(th2)), (I2 + d2^2*m2)];

     v = -K*([th1;th2;dth1;dth2] - [q1;q2;dq1;dq2]) + [ddq1;ddq2];
     tau = M*v+C*[dth1;dth2]+G;
    

% eq1 = v2*(0.4500*cos(th2) + 0.2865) - 5.3955*sin(th1) - 4.4145*sin(th1 + th2) + v1*(0.9000*cos(th2) + 1.5730) - 0.5000*dth2*(1.8000*dth1*sin(th2) + 0.9000*dth2*sin(th2))==0;
% eq2 = 0.2865*v2 - 4.4145*sin(th1 + th2) + v1*(0.4500*cos(th2) + 0.2865) + 0.4500*dth1*th1*sin(th2) == 0;   
% sol = solve([eq1, eq2],[v1]);

    dX(1) = dth1;
    dX(2) = dth2;
    dX(3) = (I2*tau(1) - I2*tau(2) + d2^2*m2*tau(1) - d2^2*m2*tau(2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) - d2*l1*m2*tau(2)*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);
    dX(4) = -(I2*tau(1) - I1*tau(2) - I2*tau(2) - d1^2*m1*tau(2) + d2^2*m2*tau(1) - d2^2*m2*tau(2) - l1^2*m2*tau(2) + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + I2*g*l1*m2*sin(th1) + d2*l1*m2*tau(1)*cos(th2) - 2*d2*l1*m2*tau(2)*cos(th2) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(m1*d1^2*d2^2*m2 + I2*m1*d1^2 - d2^2*l1^2*m2^2*cos(th2)^2 + d2^2*l1^2*m2^2 + I1*d2^2*m2 + I2*l1^2*m2 + I1*I2);

end    
