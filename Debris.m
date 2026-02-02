%% Constant
%clear;clc; close all
Earth
%% Debris
M_Debris1= 1000e3 ; % [kg]

%% orbital elements
Altitude_Debris1=36000e3; %[m] altitude at perigee 
e_Debris1= 0; %[-]  Eccentricity
i_Debris1=deg2rad(0); % [rad] Inclination
RAAN_Debris1=deg2rad(0);% [rad] Right Ascension of Ascending Nod
w_Debris1=deg2rad(0);% [rad] Argument of periapsis
v_Debris1=deg2rad(0);% [rad] Tru Anomoly

%%
a_Debris1= (Altitude_Debris1+R_Earth)/(1-e_Debris1); % [m] semi major axis
T_Debris1=2*pi*sqrt(a_Debris1^3/mu_Earth); %[S] orbirt period
tspan=linspace(0,1*T_Debris1,1000);
%%
R_Debris1= a_Debris1*(1-e_Debris1^2)/(1+e_Debris1*cos(v_Debris1));% position Debris
V_Debris1= sqrt(mu_Earth*(2/R_Debris1-1/a_Debris1)); %[m/s] orbit speed
%%
x0_Debris1= R_Debris1;
y0_Debris1=0;
z0_Debris1=0;
xdot0_Debris1=0;
ydot0_Debris1=V_Debris1*cos(i_Debris1);
zdot0_Debris1=V_Debris1*sin(i_Debris1);


Initial_states_Debris1=[x0_Debris1,y0_Debris1,z0_Debris1,xdot0_Debris1,ydot0_Debris1,zdot0_Debris1];

opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
[t_out_Debris1,state_out_Debris1] = ode45(@(t,x) Debris1(t,x,mu_Earth), tspan, Initial_states_Debris1(:), opts);

r = state_out_Debris1(:,1:3);

figure(1); hold on; grid on; axis equal
plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1.5)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
title('Debris1 two-body orbit (fixed)')

% Draw Earth
[xe,ye,ze] = sphere(50);
surf(R_Earth*xe, R_Earth*ye, R_Earth*ze, 'EdgeColor','none', 'FaceAlpha',0.2);

function data_Debris1=Debris1(~,states_Debris1,mu_Earth)

r_Debris1=states_Debris1(1:3);
v_Debris1=states_Debris1(4:6);

r_norm_Debris=norm(r_Debris1);
acc_Debris1=-mu_Earth*r_Debris1/r_norm_Debris^3;

data_Debris1=[v_Debris1;acc_Debris1];

end
