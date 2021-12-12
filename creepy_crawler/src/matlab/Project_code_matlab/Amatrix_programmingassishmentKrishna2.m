%% Krishna Madhurkar
% RBE 502 Programming Assignment 2

syms r1 r2 l1 l2 m1 m2 I1 I2 q1 q2 qdot1 qdot2 g qddot1 qddot2 F X1 X2 X3 X4 u1 u2 lambda
syms k11 k12 k13 k14 k21 k22 k23 k24

v1 = simplify((r1*qdot1*cos(q1))^2 + (-r1*qdot1*sin(q1))^2);
x2dot = (qdot1*l1*cos(q1))+(r2*(qdot1+qdot2)*cos(q1+q2));
y2dot = (-(qdot1*l1*sin(q1)+r2*(qdot1+qdot2)*sin(q1+q2)));
v2 = simplify(x2dot^2+y2dot^2);
K = simplify(0.5*m1*v1 + 0.5*m2*v2 +0.5*I1*((qdot1)^2)+ 0.5*I2*(qdot1+qdot2)^2);
P1 = m1*g*r1*cos(q1);
P2 = m2*g*(l1*cos(q1) + r2*cos(q1+q2));
P = P1+P2;
L = simplify(K-P);

q = [q1; q2];
qd = [qdot1; qdot2];
qdd = [qddot1; qddot2];

DLDq1 = jacobian(L, q1);                           
DLDqdot1 = jacobian(L, qdot1);                      

DLDq2 = simplify(jacobian(L,q2));                  
DLDqdot2 = simplify(jacobian(L,qdot2));             

dDL_dtDdotq1 = jacobian(DLDqdot1,[q;qd])*[qd;qdd];
dDL_dtDdotq2 = jacobian(DLDqdot2,[q;qd])*[qd;qdd];  

eq1 = dDL_dtDdotq1-DLDq1-u1;                           
eq2 = dDL_dtDdotq2-DLDq2-u2;                           

sol = solve([eq1==0,eq2==0],[qddot1,qddot2]);
simplify(sol.qddot1);
simplify(sol.qddot2);


%% part a 
%
X3dot = subs(eq1,{q1,q2,qdot1,qdot2,qddot1,qddot2,u1},{X1,X2,0,0,0,0,0});
X4dot = subs(eq2,{q1,q2,qdot1,qdot2,qddot1,qddot2,u2},{X1,X2,0,0,0,0,0});
[x1,x2] = solve([X3dot==0,X4dot==0],[X1,X2])

%% part b Jacob Linerization 

X3dot = subs(sol.qddot1,{q1,q2,qdot1,qdot2},{X1,X2,X3,X4});
X4dot = subs(sol.qddot2,{q1,q2,qdot1,qdot2},{X1,X2,X3,X4});
x = [X1;X2;X3;X4];
u = [u1;u2];
X_dot =[X3; X4; X3dot; X4dot];
A= jacobian(X_dot,x)
A_symb = simplify(subs(A,{m1,m2,l1,l2,r1,r2,I1,I2,g},{1,1,1,1,0.45,0.45,0.084,0.084,9.81}))

B = jacobian(X_dot,u) 
B_symb = simplify(subs(B,{m1,m2,l1,l2,r1,r2,I1,I2,g},{1,1,1,1,0.45,0.45,0.084,0.084,9.81}))

%%  linerization around equilibrium 
A = simplify(subs(A_symb,{X1,X2,X3,X4},{0,0,0,0}));
A = vpa(A,4)
B = simplify(subs(B_symb,{X1,X2,X3,X4},{0,0,0,0}));
B = vpa(B,4)

%% stability 

A1 = double(subs(A_symb,{X1,X2,X3,X4},{0,0,0,0}));

[V1,D1] = eig(A1);
vpa(D1,4) %unstable

A2 = double(subs(A_symb,{X1,X2,X3,X4},{pi,0,0,0}));

[V2,D2] = eig(A3);
vpa(D2,4) %stable at pi,0

A3 = double(subs(A_symb,{X1,X2,X3,X4},{0,pi,0,0}));

[V3,D3] = eig(A3);
vpa(D3,4) %unstable

%% part d
C=[B A*B]
matrix_rank=rank(C)  %Rank == 4 therefore full rank

%% part e

% K=[k11 k12 k13 k14;...
%    k21 k22 k23 k24];
% Acl= A-(B*K)
% Acl_pol=(det(Acl-(lambda*eye(4))))
lambda = [-1,-2,-4,-6];
A = double(A);
B = double(B);
Kn = place(A,B,lambda);
%Control Law
u= -(Kn*x)
z=(A*x)+(B*u)
[t,y] = ode45(@ode_2linkmani,[0,10],[30*pi/180;45*pi/180;0;0]);

subplot(2,2,1)
plot(t,y(:,1),'-');
xlabel('Time t');
ylabel('y1');
legend('q1(rad)');

subplot(2,2,2)
plot(t,y(:,2),'-r');
xlabel('Time t');
ylabel('y2');
legend('q2(rad)');

subplot(2,2,3)
plot(t,y(:,3),'-b');
xlabel('Time t');
ylabel('y3');
legend('qdot1(rad/s)');

subplot(2,2,4)
plot(t,y(:,4),'-g');
xlabel('Time t');
ylabel('y4');
legend('qdot2(rad/s)');
