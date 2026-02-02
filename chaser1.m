%% Constant
%clear;clc; close all
Earth
%% chaser
M_chaser1= 20 ; % [kg]

%% orbital elements
Altitude_chaser1=1000e3; %[m] altitude at perigee 
e_chaser1= 0.8; %[-]  Eccentricity
i_chaser1=deg2rad(52); % [rad] Inclination
RAAN_chaser1=deg2rad(0);% [rad] Right Ascension of Ascending Nod
w_chaser1=deg2rad(0);% [rad] Argument of periapsis
v_chaser1=deg2rad(0);% [rad] Tru Anomoly

%%
a_chaser1= (Altitude_chaser1+R_Earth)/(1-e_chaser1); % [m] semi major axis
T_chaser1=2*pi*sqrt(a_chaser1^3/mu_Earth); %[S] orbirt period
tspan=linspace(0,1*T_chaser1,1000);
%%
R_chaser1= a_chaser1*(1-e_chaser1^2)/(1+e_chaser1*cos(v_chaser1));% position chaser
V_chaser1= sqrt(mu_Earth*(2/R_chaser1-1/a_chaser1)); %[m/s] orbit speed
%%
x0_chaser1= R_chaser1;
y0_chaser1=0;
z0_chaser1=0;
xdot0_chaser1=0;
ydot0_chaser1=V_chaser1*cos(i_chaser1);
zdot0_chaser1=V_chaser1*sin(i_chaser1);


Initial_states_chaser1=[x0_chaser1,y0_chaser1,z0_chaser1,xdot0_chaser1,ydot0_chaser1,zdot0_chaser1];

opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
[t_out_chaser1,state_out_chaser1] = ode45(@(t,x) Chaser1(t,x,mu_Earth), tspan, Initial_states_chaser1(:), opts);

r = state_out_chaser1(:,1:3);

figure(1); hold on; grid on; axis equal
plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1.5)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
title('Chaser1 two-body orbit (fixed)')

% Draw Earth
[xe,ye,ze] = sphere(50);
surf(R_Earth*xe, R_Earth*ye, R_Earth*ze, 'EdgeColor','none', 'FaceAlpha',0.2);

function data_chaser1=Chaser1(~,states_chaser1,mu_Earth)

r_chaser1=states_chaser1(1:3);
v_chaser1=states_chaser1(4:6);

r_norm_chaser=norm(r_chaser1);
acc_chaser1=-mu_Earth*r_chaser1/r_norm_chaser^3;

data_chaser1=[v_chaser1;acc_chaser1];

end
