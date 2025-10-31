%INS
%% ---------------------------------------------------------------
%  Avionics Navigation Systems Visualizer
%  Author: Abdulla Sadoun
%  Course: ENGR 6461 / AERO 482
%  Ref: Rodrigues - Fundamentals of Avionics Navigation Systems (Ch.1–5)
%% ---------------------------------------------------------------

clear; clc; close all;

%% SECTION 1 — VECTOR PROJECTIONS & NORMALS  (Ch. 2.3.1)
disp('--- Vector Projection Demo ---');

% Define two arbitrary vectors in 3D
a = [3 2 1];
b = [2 0 1];

% Orthogonal projection of a onto b
proj_ab = (dot(a,b)/norm(b)^2) * b;

% Oblique projection example (projection along a direction not ⟂ plane)
k = [0 0 1]; % direction of projection (oblique)
proj_oblique = a - ((dot(a,b)/dot(b,k))*k);
fprintf('Oblique projection: [%g %g %g]\n', proj_oblique);

% Normal to the plane defined by a and b
n = cross(a,b)/norm(cross(a,b));


fprintf('Orthogonal projection: [%g %g %g]\n',proj_ab);
fprintf('Plane normal: [%g %g %g]\n',n);

% Plot (Orthogonal Projection)
figure('Name','Vector Projections');
quiver3(0,0,0,a(1),a(2),a(3),'b','LineWidth',2); hold on;
quiver3(0,0,0,b(1),b(2),b(3),'r','LineWidth',2);
quiver3(0,0,0,proj_ab(1),proj_ab(2),proj_ab(3),'g','LineWidth',2);
quiver3(0,0,0,proj_oblique(1),proj_oblique(2),proj_oblique(3),'m','LineWidth',2);
legend('Vector a','Vector b','Orthogonal projection a→b','Oblique Projection a→b direction k');
axis equal; grid on; xlabel('x'); ylabel('y'); zlabel('z');
title('Vector Projection and Normal Visualization');





%% SECTION 2 — SPHERICAL EARTH DISTANCE (Ch. 1.4)
disp('--- Great-Circle Distance Demo ---');

Re = 6371; % km
phi1 = deg2rad(44.67);  % Halifax latitude
lam1 = deg2rad(-63.58);
phi2 = deg2rad(28.61);  % New Delhi latitude
lam2 = deg2rad(77.21);

% Great-circle central angle
c = acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lam2-lam1));
gc_distance = Re * c;
fprintf('Great-circle distance: %.1f km\n', gc_distance);

% Plot great circle on a globe
figure('Name','Great Circle Visualization');
axesm('globe'); axis off; gridm on; 
plot3m([rad2deg(phi1) rad2deg(phi2)], [rad2deg(lam1) rad2deg(lam2)], [0 0], 'r','LineWidth',2);
title('Great Circle Path (Halifax–New Delhi)');
view(45,30);


%% SECTION 3 — ELLIPSOIDAL DISTANCE (optional refinement)
disp('--- Ellipsoidal Earth Distance (Vincenty approx.) ---');

a_ell = 6378.137; f = 1/298.257223563; b_ell = a_ell*(1-f);
% Use MATLAB Mapping Toolbox if available
if exist('distance','file')
    dist_deg = distance('gc', rad2deg(phi1), rad2deg(lam1), rad2deg(phi2), rad2deg(lam2));
    dist_km = deg2km(dist_deg);
    fprintf('Ellipsoidal great-circle distance (MATLAB): %.1f km\n', dist_km);
end


%% SECTION 4 — INS SENSOR ERROR VISUALIZATION (Ch. 4–5)
disp('--- INS Accelerometer / Gyro Error Propagation ---');

t = 0:0.1:100; % seconds
bias_accel = 50e-6; % 50 micro-g (in g's)
bias_gyro = 0.01 * pi/180; % 0.01 deg/s
vel_error = bias_accel * 9.80665 * t; % velocity error m/s
att_error = bias_gyro * t * (180/pi); % attitude drift deg

figure('Name','INS Error Growth');
subplot(2,1,1);
plot(t, vel_error); grid on; xlabel('Time [s]'); ylabel('Δv [m/s]');
title('Accelerometer Bias → Velocity Error Growth');
subplot(2,1,2);
plot(t, att_error); grid on; xlabel('Time [s]'); ylabel('Attitude Drift [deg]');
title('Gyro Bias → Attitude Drift');


%% SECTION 5 — ORTHOGONALITY CHECK (any 3 basis vectors)
disp('--- Orthogonality Check ---');

E = [1 0 0; 0 1 0; 0 0 1]; % ENU example (identity)
isOrthogonal = norm(E*E' - eye(3)) < 1e-10;
fprintf('Is transformation orthogonal? %d\n', isOrthogonal);

disp('Visualization complete ✅');
