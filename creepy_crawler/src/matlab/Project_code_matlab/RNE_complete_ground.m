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
g = -9.8;  % [m/s2] Gravity acceleration (aligned with the Y axis)
%syms L1 L2
%% *** STEP 1 ***
% Calculate the home configurations of each link, expressed in the space frame                
M1 = [1 0 0 L1/2; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {1} expressed in the {0} (space) reference frame
M2 = [1 0 0 L1; 0 1 0 -L2/2; 0 0 1 0; 0 0 0 1]; % pose of frame {2} expressed in the {0} (space) reference frame
M3 = [1 0 0 L1; 0 1 0 -L2; 0 0 1 0; 0 0 0 1]; % pose of frame {3} expressed in the {0} (space) reference frame

% Calculate the home configurations of each link, expressed w.r.t. the previous link frame
M01 = inv(eye(4)) * M1 ; % pose of frame {1} expressed in the {0} (space) reference frame
M12 = inv(M1) * M2; % pose of frame {2} expressed in the {1} reference frame
M23 = inv(M2) * M3; % pose of frame {3} expressed in the {2} reference frame

% Define the screw axes of each joint, expressed in the space frame
%S = zeros(6,n);
S = [0 0 1 0 0 0;
    0 0 1 -cross([0 0 1], [L1 0 0])]';
% Calculate the screw axes of each joint, expressed in the local link frame
A1 = adjoint(inv(M1)) * S(:,1)
A2 = adjoint(inv(M2)) * S(:,2)
A = [A1,A2];

% Initialize the twists and accelerations of each link
V1 = zeros(6,1);
V2 = zeros(6,1);
Vd1 = zeros(6,1);
Vd2 = zeros(6,1);
 
% Initialize the joint positions and velocities
syms q [1 2]
syms qd [1 2]
syms qdd [1 2]

%% *** STEP 2 ***
V0 = zeros(6,1);
Vd0 = [0 0 0 0 -g 0]'; 
% Forward Iteration - First Link
T01 = fkine(A(:,1), M01, q(1), 'space');
V1 = adjoint(inv(T01)) * V0 + A(:,1) *qd(1) % Link Velocity
Vd1 = adjoint(inv(T01)) * Vd0 + ad(V1) * A(:,1) * qd(1) + A(:,1) * qdd(1) % Link Acceleration
     
% Forward Iteration - Second Link
T12 = fkine(A(:,2), M12, q(2), 'space');
V2 =  adjoint(inv(T12)) * V1 + A(:,2) * qd(2) % Link Velocity
Vd2 = adjoint(inv(T12)) * Vd1 + ad(V2) * A(:,2) * qd(2) + A(:,2) * qdd(2) % Link Acceleration

%% Step 2: Initialize the Spatial Inertia Matrices

G1 = [Inertia_box(m1, 0.01, 0.01, 0.1) zeros(3,3); zeros(3,3) m1*eye(3,3)] % Spatial Inertia Matrix for Link 1
G2 = [Inertia_box(m2, 0.01, 0.01, 0.1) zeros(3,3); zeros(3,3) m2*eye(3,3)] % Spatial Inertia Matrix for Link 2

%% Step 3: Calculate the Joint Torques
syms N
N = g
F3 = [0 0 0 0 N 0]'; % Wrench applied at the end effector

% Second joint
T23 = eye(4);
F2 = adjoint(inv(T23))'*F3 + G2 * Vd2 - ad(V2)' * G2 * V2;
u2 = F2' * A(:,2);

% First joint
F1 = adjoint(inv(T12))'*F2 + G1 * Vd1 - ad(V1)' * G1 * V1;
u1 = F1' * A(:,1);

syms th1 th2 dth1 dth2 ddq ddth1 ddth2 tau1 tau2
u1 = subs(u1, [conj(qdd1), conj(qdd2), conj(qd1), conj(qd2), conj(q1), conj(q2), qdd1, qdd2, qd1, qd2, q1, q2], [ddth1, ddth2, dth1, dth2, th1, th2, ddth1, ddth2, dth1, dth2, th1, th2])
u2 = subs(u2, [conj(qdd1), conj(qdd2), conj(qd1), conj(qd2), conj(q1), conj(q2), qdd1, qdd2, qd1, qd2, q1, q2 ], [ddth1, ddth2, dth1, dth2, th1, th2, ddth1, ddth2, dth1, dth2, th1, th2])

%%______________________________________________________
% % State Space
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

% _________________________________________________________



%%______________________________________________________
% % Set Point






% _________________________________________________________

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