clear; close; clc;

% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;

while(t < 10)
    t = toc;
    % read the joint states
    jointData = receive(JointStates);
    
    % inspect the "jointData" variableinMATLAB to get familiar with itsstructure
    % design your state feedback controllerinthe following
    
    tau1.Data = - (120360503154918929209937879*X1)/19342813113834066795298816 - (15954240967328887261846177*X2)/9671406556917033397649408 - (271778644522155125246448167*X3)/38685626227668133590597632 - (8690605978639639333605539*X4)/19342813113834066795298816 ;
    tau2.Data = - (13237720256046336399771451*X1)/19342813113834066795298816 - (305437866997096155927683355*X2)/38685626227668133590597632 - (7914769523226858116048987*X3)/38685626227668133590597632 - (462268992875058155206196229*X4)/77371252455336267181195264;
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    % you can sample data here to plot at the end
end
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore

rosshutdown;