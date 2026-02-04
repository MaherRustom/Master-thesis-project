%% Constant
%clear;clc; close all
Earth
%% chaser
M_Debris= 20 ; % [kg]

%% orbital elements
Altitude_Debris=36000e3; %[m] altitude at perigee 
e_Debris= 0.8; %[-]  Eccentricity
i_Debris=deg2rad(52); % [rad] Inclination
RAAN_Debris=deg2rad(0);% [rad] Right Ascension of Ascending Nod
w_Debris=deg2rad(0);% [rad] Argument of periapsis
v_Debris=deg2rad(0);% [rad] Tru Anomoly
%% Attitude of  chaser, refrense fram Earth, Z north pole, y east, x outward
rotation_chaser=[deg2rad(54),deg2rad(21),deg2rad(18)]; %[Rad] [roll,pitch, Yaw]

q_chaser=[....
cos(rotation_chaser(1)/2)*cos(rotation_chaser(2)/2)*cos(rotation_chaser(3)/2)+ sin(rotation_chaser(1)/2)*sin(rotation_chaser(2)/2)*sin(rotation_chaser(3)/2);
sin(rotation_chaser(1)/2)*cos(rotation_chaser(2)/2)*cos(rotation_chaser(3)/2)- cos(rotation_chaser(1)/2)*sin(rotation_chaser(2)/2)*sin(rotation_chaser(3)/2);
cos(rotation_chaser(1)/2)*sin(rotation_chaser(2)/2)*cos(rotation_chaser(3)/2)+ sin(rotation_chaser(1)/2)*cos(rotation_chaser(2)/2)*sin(rotation_chaser(3)/2);
cos(rotation_chaser(1)/2)*cos(rotation_chaser(2)/2)*sin(rotation_chaser(3)/2)- sin(rotation_chaser(1)/2)*sin(rotation_chaser(2)/2)*cos(rotation_chaser(3)/2);
];

q_chaser=q_chaser/norm(q_chaser);

% wikipedia is the source: quaternion to Euler angel
rotation_chaser(1)=atan2(2*( q_chaser(1)*q_chaser(2)+q_chaser(3)*q_chaser(4)),  1-2*(q_chaser(2)^2+q_chaser(3)^2));
rotation_chaser(2)=-pi/2+2*atan2(  sqrt(1+2*( q_chaser(1)*q_chaser(3)-q_chaser(2)*q_chaser(4)))   , sqrt(1-2*( q_chaser(1)*q_chaser(3)-q_chaser(2)*q_chaser(4)))   );
rotation_chaser(3)=atan2(2*( q_chaser(1)*q_chaser(4)+q_chaser(2)*q_chaser(3)),  1-2*(q_chaser(3)^2+q_chaser(4)^2));

angular_rotation_Debris=[0,0,0]; % [rad/s] [wx,wy,wz]
%%
a_Debris= (Altitude_Debris+R_Earth)/(1-e_Debris); % [m] semi major axis
T_Debris=2*pi*sqrt(a_Debris^3/mu_Earth); %[s] orbirt period
tspan=linspace(0,1*T_Debris,1000);
%%
R_Debris= a_Debris*(1-e_Debris^2)/(1+e_Debris*cos(v_Debris));% position chaser
V_Debris= sqrt(mu_Earth*(2/R_Debris-1/a_Debris)); %[m/s] orbit speed
%%
x0_Debris= R_Debris;
y0_Debris=0;
z0_Debris=0;
xdot0_Debris=0;
ydot0_Debris=V_Debris*cos(i_Debris);
zdot0_Debris=V_Debris*sin(i_Debris);


Initial_states_Debris=[x0_Debris,y0_Debris,z0_Debris,xdot0_Debris,ydot0_Debris,zdot0_Debris];

opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
[t_out_Debris,state_out_Debris] = ode45(@(t,x) Debris(t,x,mu_Earth), tspan, Initial_states_Debris(:), opts);

r = state_out_Debris(:,1:3);

figure(1); hold on; grid on; axis equal
plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1.5)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
title('Debris two-body orbit (fixed)')

% Draw Earth
[xe,ye,ze] = sphere(50);
surf(R_Earth*xe, R_Earth*ye, R_Earth*ze, 'EdgeColor','none', 'FaceAlpha',0.2);

function data_Debris=Debris(~,states_Debris,mu_Earth)

r_Debris=states_Debris(1:3);
v_Debris=states_Debris(4:6);

r_norm_chaser=norm(r_Debris);
acc_Debris=-mu_Earth*r_Debris/r_norm_chaser^3;

data_Debris=[v_Debris;acc_Debris];

end
