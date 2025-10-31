function ins_ellipsoid_project()
% Avionics Navigation (AERO482/ENGR6461) — Deterministic strapdown INS (local nav frame)
% Task coverage:
% T1: Frames & Rz(ψ) rotation (Rodrigues §2.2); projection sanity check (P = I - nn^T, §2.3.1–2.4)
% T2: Gravity & specific force (Rodrigues §4.2.1–§4.2.5)
% T3: Strapdown mechanization (deterministic) (Rodrigues §4.3)
% T4: Bias estimation in stationary windows; drift comparison
% T5: Uncertainty & 95% error ellipse (Rodrigues §5.6)
% T6: Sanity plots & short discussion prints

%% Simulation parameters
dt   = 0.1;                 % s
Tend = 100;                 % s
t    = 0:dt:Tend;

% Truth trajectory (local nav axes: x-East, y-North, z-Up)
V   = 200;                  % m/s
R   = 5000;                 % m
omega = V/R;                % rad/s
Vz  = -3;                   % m/s descent
g   = 9.80665;              % m/s^2 (constant g accepted)

% Sensor error parameters
bg_true   = [0;0;5e-5];     % rad/s biases: roll,pitch zero; yaw = 5e-5 (given)
sigma_gyro = 1e-5;          % rad/s white noise per axis
ba_true   = [0.01;0.01;0.01]; % m/s^2 bias per axis
sigma_acc  = 0.05;            % m/s^2 white noise per axis

% Stationary windows for bias estimation (seconds)
stationary_windows = [0 10; 50 60];

% Ellipsoid sanity (not used for kinematics—local frame INS only)
wgs84 = struct('a',6378137.0,'f',1/298.257223563); % semi-major, flattening
phi_nom = deg2rad(45.5);  % nominal latitude (e.g., Montreal) for radii check
[RM, RN] = wgs84_radii(phi_nom, wgs84); %#ok<NASGU> % sanity values if needed in prints

% Assemble truth states
[truth, att] = simulate_truth(t, R, omega, Vz);

% --- T1: Frame rotation and projection sanity check ---
% Verify orthogonality Rz' * Rz = I and projection P = I - nn^T
psi_sample = att.psi(1);
Rbn = Rz(psi_sample)';           % body<-nav (here body yawed by ψ about z)
assert( norm(Rbn'*Rbn - eye(3), 'fro') < 1e-12, 'R is not orthogonal');
n = [0;0;1]; P = eye(3) - n*n.'; % Projection onto horizontal plane (§2.3.1, §2.4)

% --- Generate inertial measurements (body frame) ---
meas_no_bias   = gen_imu(truth, att, g, zeros(3,1), 0, zeros(3,1), 0, dt);           % clean
meas_with_bias = gen_imu(truth, att, g, ba_true, sigma_acc, bg_true, sigma_gyro, dt); % biased+noisy

% --- T3: Mechanization w/ and w/o bias compensation (deterministic Euler) ---
x0 = truth.r(:,1); v0 = truth.v(:,1); psi0 = att.psi(1);
% First, run naive (no bias compensation)
ins_nocomp = strapdown_ins(meas_with_bias, g, x0, v0, psi0, dt, [], []);
% Then estimate biases from stationary windows and re-run with compensation
[bg_hat, ba_hat] = estimate_biases(meas_with_bias, att, g, stationary_windows, t);
ins_comp = strapdown_ins(meas_with_bias, g, x0, v0, psi0, dt, bg_hat, ba_hat);

% --- T5: Monte Carlo for covariance and 95% error ellipse ---
MC = 200;
[xy_cov, samples_xy] = monte_carlo_cov(MC, truth, att, g, ba_true, sigma_acc, bg_true, sigma_gyro, dt, x0, v0, psi0);
[ell_axes, ell_angles] = error_ellipse(xy_cov, 0.95); % semi-axes and principal dir

% --- T6: Plots & reporting ---
make_plots(truth, ins_nocomp, ins_comp, samples_xy, xy_cov, ell_axes, ell_angles, t);

% Print a short drift discussion (signs/magnitudes)
fprintf('\n=== Drift discussion (qualitative) ===\n');
fprintf('Yaw bias bg_z = %+g rad/s -> heading integrates with constant slope, causing lateral (cross-track) drift ~ V * heading error.\n', bg_true(3));
fprintf('Accel biases ba ~= [%g %g %g] m/s^2 -> velocity error grows ~ ba * t, position ~ 0.5*ba*t^2 in each axis.\n', ba_true);
fprintf('Bias compensation estimated bg_hat = [%g %g %g], ba_hat = [%g %g %g].\n', bg_hat, ba_hat);
fprintf('RMSE (100 s): No-comp vs Comp (m):\n');
rmse_nocomp = sqrt(mean((ins_nocomp.r - truth.r).^2,2));
rmse_comp   = sqrt(mean((ins_comp.r   - truth.r).^2,2));
fprintf('  x: %.2f -> %.2f | y: %.2f -> %.2f | z: %.2f -> %.2f\n', rmse_nocomp(1), rmse_comp(1), rmse_nocomp(2), rmse_comp(2), rmse_nocomp(3), rmse_comp(3));

end

%% ========= Helpers =========

function [truth, att] = simulate_truth(t, R, omega, Vz)
% Truth trajectory and attitude (roll=pitch=0; yaw ψ=ωt)
x = R*cos(omega*t);
y = R*sin(omega*t);
z = Vz*t;

vx = -R*omega*sin(omega*t);
vy =  R*omega*cos(omega*t);
vz =  Vz + 0*t;

ax = -R*omega^2*cos(omega*t);
ay = -R*omega^2*sin(omega*t);
az = 0*t;

truth.r = [x; y; z];
truth.v = [vx; vy; vz];
truth.a = [ax; ay; az];

att.roll  = zeros(size(t));
att.pitch = zeros(size(t));
att.psi   = omega*t; % yaw increases linearly
end

function meas = gen_imu(truth, att, g, ba, sigma_acc, bg, sigma_gyro, dt)
% Generate gyro (ψdot + noise+bias) and accel specific force (body) with noise+bias
N = numel(att.psi);
% Gravity in nav (ENU, Up positive): gn = [0;0;-g]
gn = repmat([0;0;-g],1,N);

% True body rates: only yaw dot = ω, roll/pitch 0
gyro_true = [zeros(1,N); zeros(1,N); gradient(att.psi, dt)];
% Measured gyro
meas.gyro = gyro_true + bg + sigma_gyro*randn(3,N);

% Specific force in nav: fn = a_n - g_n  (Rodrigues §4.2.4)
fn = truth.a - gn;

% Rotate to body: fb = R_n^b fn
fb = zeros(3,N);
for k=1:N
    Rnb = Rz(att.psi(k)); % nav->body is yaw about z
    fb(:,k) = Rnb*fn(:,k);
end

% Measured accel
meas.acc = fb + ba + sigma_acc*randn(3,N);
meas.t   = (0:N-1)*dt;
meas.psi0 = att.psi(1);
end

function out = strapdown_ins(meas, g, x0, v0, psi0, dt, bg_hat, ba_hat)
% Deterministic Euler mechanization in local nav frame
% ψ_k = ψ_{k-1} + (ωz_meas - bg_hat_z) Δt
% f_n = R_b^n (f_b_meas - b_a) ; a_n = f_n + g_n ; v, r integrate forward
N = numel(meas.t);
psi = zeros(1,N); psi(1)=psi0;
v   = zeros(3,N); v(:,1)=v0;
r   = zeros(3,N); r(:,1)=x0;

if isempty(bg_hat), bg_hat = zeros(3,1); end
if isempty(ba_hat), ba_hat = zeros(3,1); end

for k=2:N
    % attitude
    psidot = meas.gyro(3,k-1) - bg_hat(3);
    psi(k) = psi(k-1) + psidot*dt;

    % accel: body->nav then add gravity
    fb_corr = meas.acc(:,k-1) - ba_hat;
    Rbn = Rz(psi(k))';     % body->nav
    fn  = Rbn*fb_corr;
    an  = fn + [0;0;-g];

    % integrate velocity and position
    v(:,k) = v(:,k-1) + an*dt;
    r(:,k) = r(:,k-1) + v(:,k)*dt;
end

out.psi = psi;
out.v   = v;
out.r   = r;
end

function [bg_hat, ba_hat] = estimate_biases(meas, att, g, windows, tvec)
% Estimate constant gyro & accel biases using stationary windows
% At rest (roll=pitch=0), true psidot ≈ 0; expected fb_true ≈ [0;0; +g] in body
sel = false(size(tvec));
for i=1:size(windows,1)
    sel = sel | (tvec>=windows(i,1) & tvec<=windows(i,2));
end

% Gyro bias: mean measured in windows
bg_hat = mean(meas.gyro(:,sel),2);

% Acc bias: mean(measured - expected)
g_body = zeros(3,sum(sel));
k = 1;
for j=find(sel)
    % yaw-only rotation: z-axis unaffected → expected [0;0;g]
    g_body(:,k) = [0;0;g];
    k=k+1;
end
ba_hat = mean(meas.acc(:,sel) - g_body, 2);
end

function [Kxy, samples_xy] = monte_carlo_cov(MC, truth, att, g, ba_true, sigma_acc, bg_true, sigma_gyro, dt, x0, v0, psi0)
samples_xy = zeros(2, MC);
for m=1:MC
    meas = gen_imu(truth, att, g, ba_true, sigma_acc, bg_true, sigma_gyro, dt);
    ins  = strapdown_ins(meas, g, x0, v0, psi0, dt, ba_true*0, bg_true*0); %#ok<*NASGU>
    err_xy = ins.r(1:2,end) - truth.r(1:2,end);
    samples_xy(:,m) = err_xy;
end
mu = mean(samples_xy,2);
D  = samples_xy - mu;
Kxy = (D*D.')/(MC-1);  % sample covariance
end

function [axes_len, theta] = error_ellipse(K, conf)
% Rodrigues §5.6: eigen-decompose K = W Λ W^T; 95% ellipse uses χ²_{2,0.95} ≈ 5.991
if nargin<2, conf=0.95; end
chi2 = 5.991; % 95% for 2 dof
[V, D] = eig((K+K')/2);
lambda = max(real(diag(D)),0);
axes_len = sqrt(chi2 * lambda(:));
% principal direction angle (radians from x-axis)
theta = atan2(V(2,2), V(1,2)); % angle of major axis column (choose second by convention)
end

function make_plots(truth, ins_nocomp, ins_comp, samples_xy, Kxy, ell_axes, ell_theta, t)
figure('Name','3D Truth vs INS (bias-off vs bias-comp)');
plot3(truth.r(1,:), truth.r(2,:), truth.r(3,:), 'k-', 'LineWidth',1.5); hold on;
plot3(ins_nocomp.r(1,:), ins_nocomp.r(2,:), ins_nocomp.r(3,:), 'r--');
plot3(ins_comp.r(1,:),   ins_comp.r(2,:),   ins_comp.r(3,:),   'b-');
grid on; xlabel('x_E [m]'); ylabel('y_N [m]'); zlabel('z_U [m]');
legend('Truth','INS no-comp','INS bias-comp','Location','best'); title('T6(i): 3-D trajectories');

% N–E path with 95% ellipse at t=100 s
figure('Name','NE path + 95% error ellipse at 100 s');
plot(truth.r(1,:), truth.r(2,:), 'k-', 'LineWidth',1.2); hold on;
plot(ins_nocomp.r(1,:), ins_nocomp.r(2,:), 'r--');
plot(ins_comp.r(1,:),   ins_comp.r(2,:),   'b-');
% Draw ellipse centered at truth end
cx = truth.r(1,end); cy = truth.r(2,end);
th = linspace(0,2*pi,200); a=ell_axes(2); b=ell_axes(1); % ensure a ≥ b by eigen ordering (may swap)
R = [cos(ell_theta) -sin(ell_theta); sin(ell_theta) cos(ell_theta)];
ell = R * [a*cos(th); b*sin(th)];
plot(cx+ell(1,:), cy+ell(2,:), 'm-', 'LineWidth',1.5);
axis equal; grid on; xlabel('x_E [m]'); ylabel('y_N [m]');
legend('Truth','INS no-comp','INS bias-comp','95% error ellipse','Location','best');
title('T6(ii): N–E path and 95% error ellipse at 100 s (Monte Carlo)');

% Time histories of position error components
figure('Name','Position errors vs time');
err_nc = ins_nocomp.r - truth.r;
err_c  = ins_comp.r   - truth.r;
subplot(3,1,1); plot(t,err_nc(1,:),'r--', t,err_c(1,:),'b-'); grid on; ylabel('ex [m]'); title('T6(iii): Position errors'); legend('no-comp','bias-comp');
subplot(3,1,2); plot(t,err_nc(2,:),'r--', t,err_c(2,:),'b-'); grid on; ylabel('ey [m]');
subplot(3,1,3); plot(t,err_nc(3,:),'r--', t,err_c(3,:),'b-'); grid on; ylabel('ez [m]'); xlabel('t [s]');

% Scatter of MC terminal errors
figure('Name','Monte Carlo terminal errors (x,y)');
plot(samples_xy(1,:), samples_xy(2,:), 'o'); grid on; axis equal;
xlabel('x error [m]'); ylabel('y error [m]');
title(sprintf('T5: MC terminal errors, Kxy = [%.2f %.2f; %.2f %.2f]', Kxy(1,1),Kxy(1,2),Kxy(2,1),Kxy(2,2)));
end

%% --- Math utilities (Rodrigues §2.2, §3.2.2 sanity, etc.) ---

function R = Rz(psi)
% Rotation about +z (yaw): nav->body when roll=pitch=0 and yaw=psi
c = cos(psi); s = sin(psi);
R = [ c  s  0;
     -s  c  0;
      0  0  1];
end

function [RM, RN] = wgs84_radii(phi, wgs)
% Meridional (RM) and prime-vertical (RN) radii of curvature on the ellipsoid
% (sanity values only; INS remains local)
a = wgs.a; f = wgs.f; e2 = f*(2-f);
den = (1 - e2*sin(phi).^2);
RN  = a ./ sqrt(den);
RM  = a*(1 - e2) ./ (den.^(3/2));
end