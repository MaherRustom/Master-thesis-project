%% Constant
%clear;clc; close all
Earth
%% chaser1
M_chaser11= 20 ; % [kg]

%% orbital elements
Altitude_chaser11=6000e3; %[m] altitude at perigee 
e_chaser11= 0.9; %[-]  Eccentricity
i_chaser11=deg2rad(45); % [rad] Inclination
RAAN_chaser11=deg2rad(90);% [rad] Right Ascension of Ascending Nod
w_chaser11=deg2rad(0);% [rad] Argument of periapsis
v_chaser11=deg2rad(0);% [rad] Tru Anomoly

%% Convert orbital elements to Cartesian state
r_p_chaser11 = R_Earth + Altitude_chaser11;     % perigee radius
a_chaser11   = r_p_chaser11/(1 - e_chaser11);   % semi-major axis

T_chaser11 = 2*pi*sqrt(a_chaser11^3/mu_Earth);
tspan = linspace(0, T_chaser11, 1000);

% Position and velocity in perifocal frame
r_p_chaser11=a_chaser11*(1-e_chaser11^2)/(1+e_chaser11*cos(v_chaser11));
r_pf = r_p_chaser11*...
       [cos(v_chaser11); sin(v_chaser11); 0];

v_pf = sqrt(mu_Earth*(2/r_p_chaser11-1/a_chaser11)) * ...
       [-sin(v_chaser11); cos(v_chaser11); 0];

% Rotation matrices
R3_W = [ cos(RAAN_chaser11) -sin(RAAN_chaser11) 0;
         sin(RAAN_chaser11)  cos(RAAN_chaser11) 0;
         0                   0                  1];

R1_i = [1 0 0;
        0 cos(i_chaser11) -sin(i_chaser11);
        0 sin(i_chaser11)  cos(i_chaser11)];

% R1_i = [cos(i_chaser11)  0      sin(i_chaser11);
%         0                1                 0;
%         -sin(i_chaser11) 0    cos(i_chaser11)];

R3_w = [ cos(w_chaser11) -sin(w_chaser11) 0;
         sin(w_chaser11)  cos(w_chaser11) 0;
         0                0               1];

Q_pX = R3_W * R1_i * R3_w;

% Initial state in inertial frame
r0_chaser11 = Q_pX * r_pf;
v0_chaser11 = Q_pX * v_pf;

Initial_states_chaser11 = [r0_chaser11; v0_chaser11];

opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
[t_out_chaser11,state_out_chaser11] = ode45(@(t,x) chaser11(t,x,mu_Earth), tspan, Initial_states_chaser11(:), opts);

r = state_out_chaser11(:,1:3);

figure(1); hold on; grid on; axis equal
plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1.5)
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
title('chaser11 two-body orbit (fixed)')

% Draw Earth
[xe,ye,ze] = sphere(50);
surf(R_Earth*xe, R_Earth*ye, R_Earth*ze, 'EdgeColor','none', 'FaceAlpha',0.2);

%% Attitude of  chaser1, refrense fram Earth, Z north pole, y east, x outward
rotation_chaser1=[deg2rad(0),deg2rad(0),deg2rad(0)]; %[Rad] [roll,pitch, Yaw]

q_chaser1=[....
cos(rotation_chaser1(1)/2)*cos(rotation_chaser1(2)/2)*cos(rotation_chaser1(3)/2)+ sin(rotation_chaser1(1)/2)*sin(rotation_chaser1(2)/2)*sin(rotation_chaser1(3)/2);
sin(rotation_chaser1(1)/2)*cos(rotation_chaser1(2)/2)*cos(rotation_chaser1(3)/2)- cos(rotation_chaser1(1)/2)*sin(rotation_chaser1(2)/2)*sin(rotation_chaser1(3)/2);
cos(rotation_chaser1(1)/2)*sin(rotation_chaser1(2)/2)*cos(rotation_chaser1(3)/2)+ sin(rotation_chaser1(1)/2)*cos(rotation_chaser1(2)/2)*sin(rotation_chaser1(3)/2);
cos(rotation_chaser1(1)/2)*cos(rotation_chaser1(2)/2)*sin(rotation_chaser1(3)/2)- sin(rotation_chaser1(1)/2)*sin(rotation_chaser1(2)/2)*cos(rotation_chaser1(3)/2);
];

q_chaser1=q_chaser1/norm(q_chaser1);

% wikipedia is the source: quaternion to Euler angel
rotation_chaser1(1)=atan2(2*( q_chaser1(1)*q_chaser1(2)+q_chaser1(3)*q_chaser1(4)),  1-2*(q_chaser1(2)^2+q_chaser1(3)^2));
rotation_chaser1(2)=-pi/2+2*atan2(  sqrt(1+2*( q_chaser1(1)*q_chaser1(3)-q_chaser1(2)*q_chaser1(4)))   , sqrt(1-2*( q_chaser1(1)*q_chaser1(3)-q_chaser1(2)*q_chaser1(4)))   );
rotation_chaser1(3)=atan2(2*( q_chaser1(1)*q_chaser1(4)+q_chaser1(2)*q_chaser1(3)),  1-2*(q_chaser1(3)^2+q_chaser1(4)^2));



angular_rotation_chaser11=[0,0,0]; % [rad/s] [wx,wy,wz]

%% Functions 

function data_chaser11=chaser11(~,states_chaser11,mu_Earth)

r_chaser11=states_chaser11(1:3);
v_chaser11=states_chaser11(4:6);

r_norm_chaser1=norm(r_chaser11);
acc_chaser11=-mu_Earth*r_chaser11/r_norm_chaser1^3;

data_chaser11=[v_chaser11;acc_chaser11];

end

function q = eulerToQuat(roll, pitch, yaw)
    % Converts roll, pitch, yaw to quaternion
    % Rotation sequence: ZYX
    % roll  = rotation about x-axis
    % pitch = rotation about y-axis
    % yaw   = rotation about z-axis
    %
    % Output format:
    % q = [w x y z]

    cr = cos(roll/2);
    sr = sin(roll/2);

    cp = cos(pitch/2);
    sp = sin(pitch/2);

    cy = cos(yaw/2);
    sy = sin(yaw/2);

    w = cr*cp*cy + sr*sp*sy;
    x = sr*cp*cy - cr*sp*sy;
    y = cr*sp*cy + sr*cp*sy;
    z = cr*cp*sy - sr*sp*cy;

    q = [w x y z];
    q=q/norm(q);
end


function q = quat_multiply(a,b)

    a0 = a(1); a1 = a(2); a2 = a(3); a3 = a(4);
    b0 = b(1); b1 = b(2); b2 = b(3); b3 = b(4);

    q = [
        a0*b0 - a1*b1 - a2*b2 - a3*b3;
        a0*b1 + a1*b0 + a2*b3 - a3*b2;
        a0*b2 - a1*b3 + a2*b0 + a3*b1;
        a0*b3 + a1*b2 - a2*b1 + a3*b0
    ];

end

function q_inv = quat_inverse(q)

    q_conj = [q(1); -q(2); -q(3); -q(4)];
    q_inv = q_conj / dot(q,q);

end
