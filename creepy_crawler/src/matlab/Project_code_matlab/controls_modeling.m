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
syms u1 u2 tau1 tau2

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

solution = solve([u1==tau1 u2==tau2],[ddth1 ddth2]);
dX3 = solution.ddth1
dX4 = solution.ddth2

% dX = sym('dX', [4,1]);
% dX(1) = X(3);
% dX(2) = X(4);
% dX(3) = solution.ddth1
% dX(4) = solution.ddth2
% 
% th1 = deg2rad(30);
% th2 = deg2rad(45);
% [t, y] = ode45(@myode, [0,10],[th1,th2,0,0]);
% plot(t,y(:,1));
% figure
% plot(t,y(:,2));
% figure
% plot(t,y(:,3));
% figure
% plot(t,y(:,4));
% 
% function dX = myode(t,X)
%     m1 = 1;
%     m2 = 1; 
%     l1 = 1; 
%     l2 = 1;
%     d1 = 0.45;
%     d2 = 0.45;
%     I1 = 0.084;
%     I2 = 0.084;
%     g = 9.81;
%     
%     dX = zeros(4,1);
%     X = num2cell(X);
%     [th1, th2, dth1, dth2] = deal(X{:});
% %     if abs(th1)>2*pi
% %         th1 = mod(th1,2*pi);
% %     end
% %     
% %     if abs(th2)>2*pi
% %         th2 = mod(th2,2*pi);
% %     end
%        
%     dX(1) = dth1;
%     dX(2) = dth2;
%     dX(3) = (g*m2*d2*sin(th1 + th2) + g*l1*m2*sin(th1) + g*m1*d1*sin(th1))/(m2*l1^2 + 2*m2*cos(th2)*l1*d2 + m1*d1^2 + m2*d2^2 + I1 + I2);
%     dX(4) = (- l1*m2*d2*sin(th2)*dth1^2 + g*m2*d2*sin(th1 + th2))/(m2*d2^2 + I2);
%     
% end  
    


