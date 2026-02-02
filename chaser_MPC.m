clear; clc; close all;
Earth;

%% MPC
Ts= 0.1; % [s] sample time
T= 1000; % [s] simulation time
steps= T/Ts; 
L= 10; %[s] MPC future calculation
d= 100; %[m] safty distance

x0=[0,0,0,0,0,0]; % chaser (1-3 xyz location) (4-6 xyz velocity)
x=[0,0,0,0,0,0]; % Debris (1-3 xyz location) (4-6 xyz velocity)
q= [1,0,0,0]'; % quaternions 
v=1; % [m/s]

Q = diag([1 1 1 1 1 1]);
R = 1*eye(3);

A=zeros(8,8)
B= zeros(8,4)

Ad = 0;
Bd = 0;