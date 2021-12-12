clc
clear all
close all
sympref('FloatingPointOutput',true);
syms u1 u2 u1_swing u2_swing
syms theta1_ddot theta2_ddot theta1_dot theta2_dot theta1 theta2
syms m1 m2 I1 I2 d1 d2 l1 l2 g

% u1= (168374578068624947*theta1_ddot)/28823037615171174400 - (theta1_ddot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)))/20 + ((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2))*(theta1_ddot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)) - theta1_ddot/20 + (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(1801439850948198410*(cos(theta2)^2 + sin(theta2)^2)) + 49/5) + (180143985094819840*sin(theta2)*((1125899906842624*theta1_ddot*sin(theta2))/(67553994410557440375*(cos(theta2)^2 + sin(theta2)^2)) + (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(10808639105689190460000*(cos(theta2)^2 + sin(theta2)^2))))/(180143985094819841*(cos(theta2)^2 + sin(theta2)^2)) + ((180143985094819841*cos(theta2) - sin(theta2))*((970375599710763*theta1_ddot*(180143985094819841*cos(theta2) - sin(theta2)))/(207691874341393106294141357775650816*(cos(theta2)^2 + sin(theta2)^2)) - (1801439850948198641*theta1_dot*theta2_dot*sin(theta2))/(1080863910568919046000*(cos(theta2)^2 + sin(theta2)^2))))/(180143985094819841*(cos(theta2)^2 + sin(theta2)^2)) - (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(36028797018963968200*(cos(theta2)^2 + sin(theta2)^2)) - 49/100;
% u2= (96316984030697011*theta2_ddot)/28823037615171174400 + (441352763482308608*sin(theta2))/(900719925474099205*(cos(theta2)^2 + sin(theta2)^2)) + (theta1_dot*(theta1_dot/20 - theta1_dot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)))*(180143985094819841*cos(theta2) - sin(theta2)))/(3602879701896396820*(cos(theta2)^2 + sin(theta2)^2)) + (1783425452438716657*theta1_dot^2*sin(theta2)*(180143985094819841*cos(theta2) - sin(theta2)))/(389422264390112076463242866967183372000*(cos(theta2)^2 + sin(theta2)^2)^2);
u1_swing= (168374578068624947*theta1_ddot)/28823037615171174400 - (theta1_ddot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)))/20 + ((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2))*(theta1_ddot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)) - theta1_ddot/20 + (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(1801439850948198410*(cos(theta2)^2 + sin(theta2)^2)) + 49/5) + (180143985094819840*sin(theta2)*((1125899906842624*theta1_ddot*sin(theta2))/(67553994410557440375*(cos(theta2)^2 + sin(theta2)^2)) + (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(10808639105689190460000*(cos(theta2)^2 + sin(theta2)^2))))/(180143985094819841*(cos(theta2)^2 + sin(theta2)^2)) + ((180143985094819841*cos(theta2) - sin(theta2))*((970375599710763*theta1_ddot*(180143985094819841*cos(theta2) - sin(theta2)))/(207691874341393106294141357775650816*(cos(theta2)^2 + sin(theta2)^2)) - (1801439850948198641*theta1_dot*theta2_dot*sin(theta2))/(1080863910568919046000*(cos(theta2)^2 + sin(theta2)^2))))/(180143985094819841*(cos(theta2)^2 + sin(theta2)^2)) - (theta1_dot*theta2_dot*(180143985094819841*cos(theta2) - sin(theta2)))/(36028797018963968200*(cos(theta2)^2 + sin(theta2)^2)) - 49/100;
u2_swing= (96316984030697011*theta2_ddot)/28823037615171174400 + (441352763482308608*sin(theta2))/(900719925474099205*(cos(theta2)^2 + sin(theta2)^2)) + (theta1_dot*(theta1_dot/20 - theta1_dot*((9007199254740992*sin(theta2)*(sin(theta2) - 180143985094819841*cos(theta2) + 360287970189639681*cos(theta2)^2 + 360287970189639681*sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2) - (9007199254740992*(180143985094819841*cos(theta2) - sin(theta2))*(cos(theta2)^2 - sin(theta2) + sin(theta2)^2))/(32451855365842673038603572247265281*(cos(theta2)^2 + sin(theta2)^2)^2)))*(180143985094819841*cos(theta2) - sin(theta2)))/(3602879701896396820*(cos(theta2)^2 + sin(theta2)^2)) + (1783425452438716657*theta1_dot^2*sin(theta2)*(180143985094819841*cos(theta2) - sin(theta2)))/(389422264390112076463242866967183372000*(cos(theta2)^2 + sin(theta2)^2)^2);

X = sym('X', [4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot;
X(4) = theta2_dot;

solution = solve([u1_swing==0 u2_swing ==0],[theta1_ddot theta2_ddot]);
disp(simplify(solution.theta1_ddot))
disp(simplify(solution.theta2_ddot))

% solution = solve([u1==0 u2==0],[theta1_ddot theta2_ddot]);
% disp(simplify(solution.theta1_ddot))
% disp(simplify(solution.theta2_ddot))

dX = sym('dX', [4,1]);
dX(1) = X(3);
dX(2) = X(4);
dX(3) = solution.theta1_ddot;
dX(4) = solution.theta2_ddot;

q1 = deg2rad(30);
q2 = deg2rad(45);
[t, y] = ode45(@myode_swing, [0,10],[q1,q2,0,0]);
plot(t,y(:,1));
figure
plot(t,y(:,2));
figure
plot(t,y(:,3));
figure
plot(t,y(:,4));

function dX = myode(t,X)
   
    dX = zeros(4,1);
    X = num2cell(X);
    [theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});
    
    if abs(theta1)>2*pi
        theta1 = mod(theta1,2*pi);
    end
    
    if abs(theta2)>2*pi
        theta2 = mod(theta2,2*pi);
    end
       
    dX(1) = theta1_dot;
    dX(2) = theta2_dot;
    dX(3) = -(5.7646e+15*(4.7704e+37*cos(theta2) - 9.5408e+37*sin(theta2) + 4.5126e+18*theta1_dot*theta2_dot + 4.8678e+35*theta1_dot*theta2_dot*cos(theta2) - 2.7022e+18*theta1_dot*theta2_dot*sin(theta2) + 4.8678e+35*theta1_dot*theta2_dot*cos(theta2)^2 - 8.1292e+35*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + 4.7704e+37))/(5.6122e+51*sin(theta2) - 2.8061e+51*cos(theta2) + 5.6122e+51*cos(theta2)*sin(theta2) + 3.7461e+51*cos(theta2)^2 - 8.9000e+51);
    dX(4) = 4.1530e-18*theta1_dot^2*sin(theta2) - 0.7481*theta1_dot^2*cos(theta2) - 146.6334*sin(theta2) - 0.7481*theta1_dot^2*cos(theta2)^2 - 6.9354e-18*theta1_dot^2 + 1.2494*theta1_dot^2*cos(theta2)*sin(theta2);

end  

function dX = myode_swing(t,X)
   
    dX = zeros(4,1);
    X = num2cell(X);
    [theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});
    
    if abs(theta1)>2*pi
        theta1 = mod(theta1,2*pi);
    end
    
    if abs(theta2)>2*pi
        theta2 = mod(theta2,2*pi);
    end
       
    dX(1) = theta1_dot;
    dX(2) = theta2_dot;
    dX(3) = -(5.7646e+15*(4.7704e+37*cos(theta2) - 9.5408e+37*sin(theta2) + 4.5126e+18*theta1_dot*theta2_dot + 4.8678e+35*theta1_dot*theta2_dot*cos(theta2) - 2.7022e+18*theta1_dot*theta2_dot*sin(theta2) + 4.8678e+35*theta1_dot*theta2_dot*cos(theta2)^2 - 8.1292e+35*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + 4.7704e+37))/(5.6122e+51*sin(theta2) - 2.8061e+51*cos(theta2) + 5.6122e+51*cos(theta2)*sin(theta2) + 3.7461e+51*cos(theta2)^2 - 8.9000e+51);
    dX(4) = 4.1530e-18*theta1_dot^2*sin(theta2) - 0.7481*theta1_dot^2*cos(theta2) - 146.6334*sin(theta2) - 0.7481*theta1_dot^2*cos(theta2)^2 - 6.9354e-18*theta1_dot^2 + 1.2494*theta1_dot^2*cos(theta2)*sin(theta2);

end  