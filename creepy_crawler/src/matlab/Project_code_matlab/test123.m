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
% syms L1 L2
%% *** STEP 1 ***
% Calculate the home configurations of each link, expressed in the space frame                
M1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {1} expressed in the {0} (space) reference frame
M2 = [1 0 0 L1; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {2} expressed in the {0} (space) reference frame
M3 = [1 0 0 L1+L2; 0 1 0 0; 0 0 1 0; 0 0 0 1]; % pose of frame {3} expressed in the {0} (space) reference frame

% Calculate the home configurations of each link, expressed w.r.t. the previous link frame
M01 = pinv(eye(4)) * M1  % pose of frame {1} expressed in the {0} (space) reference frame
M12 = pinv(M1) * M2; % pose of frame {2} expressed in the {1} reference frame
M23 = pinv(M2) * M3 % pose of frame {3} expressed in the {2} reference frame

% Define the screw axes of each joint, expressed in the space frame
%S = zeros(6,n);
S = [0 0 1 0 0 0;
    0 0 1 0 -L1 0]';

% Calculate the screw axes of each joint, expressed in the local link frame
A1 = adjoint(inv(M1)) * S(:,1)
A2 = adjoint(inv(M2)) * S(:,2)
A = [A1,A2]


function AdT = adjoint(T)
    
    R = T(1:3,1:3)
    p = T(1:3,4)
    skew(p)
    AdT = [R zeros(3); skew(p)*R R];
end