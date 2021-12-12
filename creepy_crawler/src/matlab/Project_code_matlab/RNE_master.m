clc
clear all
close all
sympref('FloatingPointOutput',true);
%% Step 1: COPY-PASTE YOUR SOLUTION FOR THE JOINT VELOCITIES AND ACCELERATIONS BELOW
%% Robot Definition:
n = 2;    % Number of links in the kinematic chain
L1 = 0.1; % [m] Length of the first link
L2 = 0.1; % [m] Length of the second link
m1 = 1;   % [kg] Mass of the first link
m2 = 1;   % [kg] Mass of the second link
g = -9.8;
% g = 9.8;  % [m/s2] Gravity acceleration (aligned with the Y axis)
% syms L1 L2
%% *** STEP 1 ***
% Calculate the home configurations of each link, expressed in the space frame                
M1 = [1 0 0 L1/2; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {1} expressed in the {0} (space) reference frame
M2 = [1 0 0 L1; 0 1 0 0; 0 0 1 -L2/2; 0 0 0 1]; % pose of frame {2} expressed in the {0} (space) reference frame
M3 = [1 0 0 L1; 0 1 0 0; 0 0 1 -L2; 0 0 0 1]; % pose of frame {3} expressed in the {0} (space) reference frame

% Calculate the home configurations of each link, expressed w.r.t. the previous link frame
M01 = pinv(eye(4)) * M1 ; % pose of frame {1} expressed in the {0} (space) reference frame
M12 = pinv(M1) * M2; % pose of frame {2} expressed in the {1} reference frame
M23 = pinv(M2) * M3; % pose of frame {3} expressed in the {2} reference frame

% Define the screw axes of each joint, expressed in the space frame
%S = zeros(6,n);
S = [0 0 1 0 0 0;
    0 1 0 -cross([0 1 0], [L1 0 0])]';
% Calculate the screw axes of each joint, expressed in the local link frame
A1 = adjoint(inv(M1)) * S(:,1);
A2 = adjoint(inv(M2)) * S(:,2);
A = [A1,A2];

% Initialize the twists and accelerations of each link
V1 = zeros(6,1);
V2 = zeros(6,1);
Vd1 = zeros(6,1);
Vd2 = zeros(6,1);
 
% Initialize the joint positions and velocities
syms q [1 2]
syms dq [1 2]
syms ddq [1 2]

%% *** STEP 2 ***
V0 = zeros(6,1);
Vd0 = [0 0 0 0 0 -g].'; 
% Forward Iteration - First Link
T01 = fkine(A(:,1), M01, q(1), 'space');
V1 = adjoint(inv(T01)) * V0 + A(:,1) *dq(1); % Link Velocity
Vd1 = adjoint(inv(T01)) * Vd0 + ad(V1) * A(:,1) * dq(1) + A(:,1) * ddq(1); % Link Acceleration
     
% Forward Iteration - Second Link
T12 = fkine(A(:,2), M12, q(2), 'space');
V2 =  adjoint(inv(T12)) * V1 + A(:,2) * dq(2); % Link Velocity
Vd2 = adjoint(inv(T12)) * Vd1 + ad(V2) * A(:,2) * dq(2) + A(:,2) * ddq(2); % Link Acceleration

%% Step 2: Initialize the Spatial Inertia Matrices

G1 = [Inertia_box(m1, 0.01, 0.01, 0.1) zeros(3,3); zeros(3,3) m1*eye(3,3)]; % Spatial Inertia Matrix for Link 1
G2 = [Inertia_box(m2, 0.01, 0.01, 0.1) zeros(3,3); zeros(3,3) m2*eye(3,3)]; % Spatial Inertia Matrix for Link 2

%% Step 3: Calculate the Joint Torques

F3_ground = [0 0 0 0 g 0]'; % Wrench applied at the end effector
F3_swing = [0 0 0 0 0 0]';
%% Calculating joint torques for ground
% Second joint
T23 = eye(4);
F2 = adjoint(inv(T23))'*F3_ground + G2 * Vd2 - ad(V2)' * G2 * V2;
u2 = F2' * A(:,2);

% First joint
F1 = adjoint(inv(T12))'*F2 + G1 * Vd1 - ad(V1)' * G1 * V1;
u1 = F1' * A(:,1);

%% Calculating joint torques for swing
T23_swing = eye(4);
F2_swing = adjoint(inv(T23_swing))'*F3_ground + G2 * Vd2 - ad(V2)' * G2 * V2;
u2_swing = F2_swing' * A(:,2);

% First joint
F1_swing = adjoint(inv(T12))'*F2 + G1 * Vd1 - ad(V1)' * G1 * V1;
u1_swing = F1_swing' * A(:,1);


u1 = subs(u1, [conj(ddq1), conj(ddq2), conj(dq1), conj(dq2), conj(q1), conj(q2), ddq1, ddq2, dq1, dq2, q1, q2], [ddq1, ddq2, dq1, dq2, q1, q2, ddq1, ddq2, dq1, dq2, q1, q2]);
u2 = subs(u2, [conj(ddq1), conj(ddq2), conj(dq1), conj(dq2), conj(q1), conj(q2), ddq1, ddq2, dq1, dq2, q1, q2], [ddq1, ddq2, dq1, dq2, q1, q2, ddq1, ddq2, dq1, dq2, q1, q2]);

u1_swing = subs(u1_swing, [conj(ddq1), conj(ddq2), conj(dq1), conj(dq2), conj(q1), conj(q2), ddq1, ddq2, dq1, dq2, q1, q2], [ddq1, ddq2, dq1, dq2, q1, q2, ddq1, ddq2, dq1, dq2, q1, q2]);
u2_swing = subs(u2_swing, [conj(ddq1), conj(ddq2), conj(dq1), conj(dq2), conj(q1), conj(q2), ddq1, ddq2, dq1, dq2, q1, q2], [ddq1, ddq2, dq1, dq2, q1, q2, ddq1, ddq2, dq1, dq2, q1, q2]);

fprintf("u1: %s\n", u1)
fprintf("u2: %s\n", u2)

% q1 = th1;
% q2 = th2;
% dq1 = dth1; 
% dq2 = dth2; 
% ddq1 = ddth1; 
% ddq2 = ddth2;

X = sym('X', [4,1]);
X(1) = q1;
X(2) = q2;
X(3) = dq1;
X(4) = dq2;

solution = solve([u1==0 u2 ==0],[ddq1 ddq2]);
disp(simplify(solution.ddq1))
disp(simplify(solution.ddq2))

dX = sym('dX', [4,1]);
dX(1) = X(3);
dX(2) = X(4);
dX(3) = solution.ddq1;
dX(4) = solution.ddq2;

q1 = deg2rad(30);
q2 = deg2rad(45);
[t, y] = ode45(@myode, [0,10],[q1,q2,0,0]);
plot(t,y(:,1));
figure
plot(t,y(:,2));
figure
plot(t,y(:,3));
figure
plot(t,y(:,4));

function dX = myode(t,X)
%     m1 = 1;
%     m2 = 1; 
%     l1 = 1; 
%     l2 = 1;
%     r1 = 0.45;
%     r2 = 0.45;
%     I1 = 0.084;
%     I2 = 0.084;
%     g = 9.81;
    
    dX = zeros(4,1);
    X = num2cell(X);
    [q1, q2, dq1, dq2] = deal(X{:});
    
    if abs(q1)>2*pi
        q1 = mod(q1,2*pi);
    end
    
    if abs(q2)>2*pi
        q2 = mod(q2,2*pi);
    end
       
    dX(1) = dq1;
    dX(2) = dq2;
    dX(3) = (1.4412e+16*(9.5408e+37*cos(q2)^3 + 9.5408e+37*cos(q2)^4 - 1.9082e+38*sin(q2)^3 + 9.5408e+37*sin(q2)^4 + 9.5408e+37*cos(q2)*sin(q2)^2 - 1.9082e+38*cos(q2)^2*sin(q2) + 1.9082e+38*cos(q2)^2*sin(q2)^2 + 9.7356e+35*dq1*dq2*cos(q2)^2 + 9.7356e+35*dq1*dq2*cos(q2)^3 + 9.0252e+18*dq1*dq2*sin(q2)^2 - 5.4043e+18*dq1*dq2*sin(q2)^3 + 9.7356e+35*dq1*dq2*cos(q2)*sin(q2)^2 - 5.4043e+18*dq1*dq2*cos(q2)^2*sin(q2) - 1.6258e+36*dq1*dq2*cos(q2)*sin(q2)))/(9.3770e+51*cos(q2)^2 - 2.8061e+52*cos(q2)*sin(q2) + 1.4030e+52*cos(q2)^3 + 1.6392e+52*cos(q2)^4 + 2.8108e+52*sin(q2)^2 - 2.8061e+52*sin(q2)^3 + 1.6392e+52*sin(q2)^4 + 1.4030e+52*cos(q2)*sin(q2)^2 - 2.8061e+52*cos(q2)^2*sin(q2) + 3.2784e+52*cos(q2)^2*sin(q2)^2);
    dX(4) = -(7.6845e-37*(1.9082e+38*sin(q2)^3 + 9.7356e+35*dq1^2*cos(q2)^2 + 9.7356e+35*dq1^2*cos(q2)^3 + 9.0252e+18*dq1^2*sin(q2)^2 - 5.4043e+18*dq1^2*sin(q2)^3 + 1.9082e+38*cos(q2)^2*sin(q2) + 9.7356e+35*dq1^2*cos(q2)*sin(q2)^2 - 5.4043e+18*dq1^2*cos(q2)^2*sin(q2) - 1.6258e+36*dq1^2*cos(q2)*sin(q2)))/(cos(q2)^2 + sin(q2)^2)^2;
        
end  
    

function AdT = adjoint(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    AdT = [R zeros(3); skew(p)*R R];
end

function I = Inertia_box(m,h,w,l)
    Ixx = m * (w^2 + h^2)/12;
    Iyy = m * (l^2 + h^2)/12;
    Izz = m * (w^2 + l^2)/12;
    I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
end