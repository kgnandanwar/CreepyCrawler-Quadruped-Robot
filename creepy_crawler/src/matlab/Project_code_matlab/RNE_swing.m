clc
clear all
close all
%% Step 1: COPY-PASTE YOUR SOLUTION FOR THE JOINT VELOCITIES AND ACCELERATIONS BELOW
%% Robot Definition:
n = 2;    % Number of links in the kinematic chain
L1 = 0.1; % [m] Length of the first link
L2 = 0.1; % [m] Length of the second link
m1 = 1;   % [kg] Mass of the first link
m2 = 1;   % [kg] Mass of the second link
g = -9.8;  % [m/s2] Gravity acceleration (aligned with the Y axis)
sympref('FloatingPointOutput', true);
% syms L1 L2
%% *** STEP 1 ***
% Calculate the home configurations of each link, expressed in the space frame                
M1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {1} expressed in the {0} (space) reference frame
M2 = [1 0 0 L1; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {2} expressed in the {0} (space) reference frame
M3 = [1 0 0 L1+L2; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {3} expressed in the {0} (space) reference frame

% Calculate the home configurations of each link, expressed w.r.t. the previous link frame
M01 = inv(eye(4)) * M1 ; % pose of frame {1} expressed in the {0} (space) reference frame
M12 = inv(M1) * M2; % pose of frame {2} expressed in the {1} reference frame
M23 = inv(M2) * M3; % pose of frame {3} expressed in the {2} reference frame

% Define the screw axes of each joint, expressed in the space frame
%S = zeros(6,n);
S = [0 0 1 0 0 0;
    0 0 1 0 -L1 0]';
% Calculate the screw axes of each joint, expressed in the local link frame
A1 = adjoint(inv(M1)) * S(:,1)
A2 = adjoint(inv(M2)) * S(:,2)
A = [A1,A2]
% 
% Initialize the twists and accelerations of each link
V1 = zeros(6,1);
V2 = zeros(6,1);
Vd1 = zeros(6,1);
Vd2 = zeros(6,1);
 
% Initialize the joint positions and velocities
syms th [1 2]
syms dth [1 2]
syms ddth [1 2]

%% *** STEP 2 ***
V0 = zeros(6,1);
Vd0 = [0 0 0 0 -g 0]'; 
% Forward Iteration - First Link
T01 = fkine(A(:,1), M01, th(1), 'space');
V1 = adjoint(inv(T01)) * V0 + A(:,1) *dth(1) % Link Velocity
Vd1 = adjoint(inv(T01)) * Vd0 + ad(V1) * A(:,1) * dth(1) + A(:,1) * ddth(1) % Link Acceleration
     
% Forward Iteration - Second Link
T12 = fkine(A(:,2), M12, th(2), 'space');
V2 =  adjoint(inv(T12)) * V1 + A(:,2) * dth(2) % Link Velocity
Vd2 = adjoint(inv(T12)) * Vd1 + ad(V2) * A(:,2) * dth(2) + A(:,2) * ddth(2) % Link Acceleration

%% Step 2: Initialize the Spatial Inertia Matrices
G1 = [zeros(3,3) zeros(3,3); zeros(3,3) m1*eye(3,3)]; % Spatial Inertia Matrix for Link 1
G2 = [zeros(3,3) zeros(3,3); zeros(3,3) m2*eye(3,3)]; % Spatial Inertia Matrix for Link 2

%% Step 3: Calculate the Joint Torques
F3 = zeros(6,1); % Wrench applied at the end effector

% Second joint
T23 = eye(4);
F2 = adjoint(inv(T23))'*F3 + G2 * Vd2 - ad(V2)' * G2 * V2;
u2 = simplify(F2' * A(:,2))

% First joint
F1 = adjoint(inv(T12))'*F2 + G1 * Vd1 - ad(V1)' * G1 * V1;
u1 = simplify(F1' * A(:,1))

syms P P1 P2
syms vx vy v_sq
syms K K1 K2
syms L
syms m1 m2 I1 I2 d1 d2 l1 l2 q th1 th2 dq dth1 dth2 ddq ddth1 ddth2
syms g real
syms dL_dth1 dL_dth2 DL_ddth1 DL_ddth2 dDL_dtddth2 dDL_dtddth1
syms u1 u2 tau1 tau2

X = sym('X', [4,1]);
X(1) = th(1);
X(2) = th(2);
X(3) = dth(1);
X(4) = dth(2);

solution = solve([u1==tau1 u2==tau2],[ddth1 ddth2]);
dX3 = solution.ddth1
dX4 = solution.ddth2

dX = sym('dX', [4,1]);
dX(1) = X(3);
dX(2) = X(4);
dX(3) = solution.ddth1
dX(4) = solution.ddth2


function AdT = adjoint(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    AdT = [R zeros(3); skew(p)*R R];
end
