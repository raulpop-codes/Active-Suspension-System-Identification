%%
% Nume si prenume: POP Raul
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 5;
n = 10;

%% Process data (fixed, do not modify)
assert(isscalar(m) && m == fix(m) && any(m == 3:6), ...
    'm must be one of the group numbers 3, 4, 5, or 6.')
assert(isscalar(n) && n == fix(n) && n >= 1 && n <= 17, ...
    'n must be an integer between 1 and 17.')

student_id = n + 17*(m-3);
rng(32452843 + 130363*student_id, 'twister')
r = rand(11,1);

p = [
    4.00 + 2.50*r(1);
    0.17 + 0.22*r(2);
    0.055 + 0.030*r(3);
    0.080 + 0.030*r(4);
    0.80 + 0.40*r(5);
    0.80 + 0.40*r(6)
];

a1 = [ ...
    2*p(1)*p(2), ...
    p(1)*p(5)/p(6)
];

a2 = [ ...
    p(1)*p(6)/p(5), ...
    p(1)*p(5)*p(3)
];

b0 = p(4)/p(6);

u_op_region = -350 + 700*r(7);

x0_slx = [
    2*p(5)/(p(1)*p(4))*(2*(r(8) >= 0.5)-1)* ...
        (6 + 24*abs(r(8)-0.5));
    p(6)/p(4)*(p(3)*p(4)*u_op_region + (-10 + 20*r(9)))
];

%% Experiment setup (fixed, do not modify)
Ts = 2e-3;  % [s], fundamental step size

% Experiment duration
q = roots([1, a1(1), a1(2)*a2(1)]);
Tfin = ceil((16/min(abs(real(q))))/Ts)*Ts*8;

du_max = 400;
umin = u_op_region - du_max; % [N], input saturation
umax = u_op_region + du_max;

q_in = 6000/pow2(10);   % resolution on the range [-3000,3000] N

ymin = -35;    % [mm], output saturation
ymax =  35;

sigma_u = 4 + 6*r(10);        % [N]
sigma_y = 0.05 + 0.08*r(11);  % [mm]

whtn_Ts_in   = 3*Ts;
whtn_pow_in  = sigma_u^2*whtn_Ts_in/2;
whtn_seed_in = 41000 + 2*student_id;

whtn_Ts_out   = 3*Ts;
whtn_pow_out  = sigma_y^2*whtn_Ts_out;
whtn_seed_out = 41001 + 2*student_id;
q_out = (ymax-ymin)/pow2(10);

clear r p q sigma_u sigma_y student_id

%% Input setup (can be changed/replaced/deleted)
fmin = 0.1/5;
fmax = 10/5;
Ain  = 385;

%% Data acquisition (use t, u, y to perform system identification)
out = sim("suspensie_activa_chirp.slx");

t = out.tout;
u = out.u;
y = out.y;

% Saturation check
tol_u = q_in/2;
tol_y = q_out/2;
sat_u = (u <= umin + tol_u) | (u >= umax - tol_u);
sat_y = (y <= ymin + tol_y) | (y >= ymax - tol_y);
if any(sat_u(:))
    warning('Input u is saturated. Reduce the input amplitude.')
end
if any(sat_y(:))
    warning('Output y is saturated. Reduce the input amplitude.')
end

% Plot data
figure
%
ax1 = subplot(2,1,1);
plot(t,u)
xlabel('t [s]')
ylabel('u [N]')
grid
%
ax2 = subplot(2,1,2);
plot(t,y)
xlabel('t [s]')
ylabel('y [mm]')
grid
%
linkaxes([ax1,ax2],'x')

%% System identification
% Cur   s 5 - Modalitatea de calibrare a intrarii
% Curs 4 - Metodologie

i1 = 6083; % first clear max for y
i2 = 6178; % first clear max for u

K = y(i1) / u(i2);

i3 = 25421; % max rez y
i4 = 25782; % min rez y

i5 = 23800; % max u
i6 = 24187; % min u

Amply = (y(i3) - y(i4))/2;
Amplu = (u(i5) - u(i6))/2;
% 3 pct de maxim in jurul rezonantei => wk-1 < wr < wk

i7 = 24698; % first peak
i8 = 25383; % rez peak
i9 = 26140; % second peak

ToscPrev  = t(i8) - t(i7);
ToscAfter = t(i9) - t(i8);
Tosc = (ToscPrev + ToscAfter)/2;

wosc = 2*pi/Tosc;% =~ wrez
wr = wosc;
Mr = Amply / Amplu;


zeta = sqrt((1 - sqrt(1 - (K/Mr)^2) )/2);
wn = wr/ sqrt(1 - 2*zeta^2);

%% Transform
H = tf(K*wn^2, [1 2*zeta*wn wn^2]);
y_sim_tf = lsim(H,u,t);

%% State Space
A = [ 0 1;
    -wn.^2 -2*zeta*wn];
B = [0;
    K*wn.^2];
C = [1 0];
D = 0;

sys = ss(A,B,C,D);
y_sim_ss = lsim(sys,u,t, [y(1) 25]);

%% Filtering
idx_val = (i9+1):(length(t)- 1000); % Unused values

% Moving average filter
windowSize = 55;
y_mov = movmean(y,windowSize);

% Butterworth filter
fs = 1/Ts;
fc = 3;
ord = 2;
[b_filt, a_filt] = butter(ord, fc/(fs/2), 'low');
y_butter = filtfilt(b_filt, a_filt, y);
%% Performance indexes
% J = ||eps||_2 / sqrt(N)
% eps_MPN = ||eps||_2 / ||y - mean(y)||_2

% State Space - nefiltrat
J_ss_full    = norm(y_sim_ss - y) / sqrt(length(y));
eMPN_ss_full = norm(y_sim_ss - y) / norm(y - mean(y));
J_ss_val     = norm(y_sim_ss(idx_val) - y(idx_val)) / sqrt(length(idx_val));
eMPN_ss_val  = norm(y_sim_ss(idx_val) - y(idx_val)) / norm(y(idx_val) - mean(y(idx_val)));

% State Space - movmean
J_ss_mov_full    = norm(y_sim_ss - y_mov) / sqrt(length(y_mov));
eMPN_ss_mov_full = norm(y_sim_ss - y_mov) / norm(y_mov - mean(y_mov));
J_ss_mov_val     = norm(y_sim_ss(idx_val) - y_mov(idx_val)) / sqrt(length(idx_val));
eMPN_ss_mov_val  = norm(y_sim_ss(idx_val) - y_mov(idx_val)) / norm(y_mov(idx_val) - mean(y_mov(idx_val)));

% State Space - Butterworth
J_ss_butter_full    = norm(y_sim_ss - y_butter) / sqrt(length(y_butter));
eMPN_ss_butter_full = norm(y_sim_ss - y_butter) / norm(y_butter - mean(y_butter));
J_ss_butter_val     = norm(y_sim_ss(idx_val) - y_butter(idx_val)) / sqrt(length(idx_val));
eMPN_ss_butter_val  = norm(y_sim_ss(idx_val) - y_butter(idx_val)) / norm(y_butter(idx_val) - mean(y_butter(idx_val)));

% Transfer Function - nefiltrat
J_tf_full    = norm(y_sim_tf - y) / sqrt(length(y));
eMPN_tf_full = norm(y_sim_tf - y) / norm(y - mean(y));
J_tf_val     = norm(y_sim_tf(idx_val) - y(idx_val)) / sqrt(length(idx_val));
eMPN_tf_val  = norm(y_sim_tf(idx_val) - y(idx_val)) / norm(y(idx_val) - mean(y(idx_val)));

% Transfer Function - movmean
J_tf_mov_full    = norm(y_sim_tf - y_mov) / sqrt(length(y_mov));
eMPN_tf_mov_full = norm(y_sim_tf - y_mov) / norm(y_mov - mean(y_mov));
J_tf_mov_val     = norm(y_sim_tf(idx_val) - y_mov(idx_val)) / sqrt(length(idx_val));
eMPN_tf_mov_val  = norm(y_sim_tf(idx_val) - y_mov(idx_val)) / norm(y_mov(idx_val) - mean(y_mov(idx_val)));

% Transfer Function - Butterworth
J_tf_butter_full    = norm(y_sim_tf - y_butter) / sqrt(length(y_butter));
eMPN_tf_butter_full = norm(y_sim_tf - y_butter) / norm(y_butter - mean(y_butter));
J_tf_butter_val     = norm(y_sim_tf(idx_val) - y_butter(idx_val)) / sqrt(length(idx_val));
eMPN_tf_butter_val  = norm(y_sim_tf(idx_val) - y_butter(idx_val)) / norm(y_butter(idx_val) - mean(y_butter(idx_val)));

%% Plots - State Space model
figure
subplot(3,1,1)
plot(t,y_sim_ss,'r', t, y, 'g');
title('State Space model - Unfiltered')
legend('Simulated', sprintf('Measured (J=%.4f, eps_{MPN}=%.2f%%)', J_ss_full, 100*eMPN_ss_full))

subplot(3,1,2)
plot(t,y_sim_ss,'r', t, y_mov, 'g');
title('State Space model - MovMean filtered')
legend('Simulated', sprintf('Filtered (J=%.4f, eps_{MPN}=%.2f%%)', J_ss_mov_full, 100*eMPN_ss_mov_full))

subplot(3,1,3)
plot(t,y_sim_ss,'r', t, y_butter, 'g');
title('State Space model - Butterworth filtered')
legend('Simulated', sprintf('Filtered (J=%.4f, eps_{MPN}=%.2f%%)', J_ss_butter_full, 100*eMPN_ss_butter_full))

%% Plots - Transfer Function model
figure
subplot(3,1,1)
plot(t,y_sim_tf,'r', t, y, 'g');
title('Transfer Function model - Unfiltered')
legend('Simulated', sprintf('Measured (J=%.4f, eps_{MPN}=%.2f%%)', J_tf_full, 100*eMPN_tf_full))

subplot(3,1,2)
plot(t,y_sim_tf,'r', t, y_mov, 'g');
title('Transfer Function model - MovMean filtered')
legend('Simulated', sprintf('Filtered (J=%.4f, eps_{MPN}=%.2f%%)', J_tf_mov_full, 100*eMPN_tf_mov_full))

subplot(3,1,3)
plot(t,y_sim_tf,'r', t, y_butter, 'g');
title('Transfer Function model - Butterworth filtered')
legend('Simulated', sprintf('Filtered (J=%.4f, eps_{MPN}=%.2f%%)', J_tf_butter_full, 100*eMPN_tf_butter_full))
