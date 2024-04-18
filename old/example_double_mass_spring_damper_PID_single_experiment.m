%% Descrpition
% This script is ment to demonstrate a case where performing VRFT on a
% Hammerstein system may result in instability due to underestimation of
% the real gain of the system.

% Questions:
% - Why are some results negative once in a while?

%% Some general settings
clc; % close;
rng(2);

%% Parameters
Ts = 0.1;
T = 0:Ts:2500;
noise_var = 0.025;

%% Nonlinearity
alpha = 0.5;
sat = @(u) min(max(u,-alpha),alpha);
f = @(u) sat(u);

%% System
% Transfer function open loop
m1 = 1; m2 = 0.5;
c1 = 0.2; c2 = 0.5;
k1 = 1; k2 = 0.5;

s = tf('s');
G = minreal(1/((m1*s^2+(c1+c2)*s+k1+k2)*(m2*s^2+c2*s+k2)-(k2+c2*s)^2));

% State space representation (augmented to include PID)
sys = ss(G);
A = sys.A; B = sys.B; C = sys.C;

% Funcitons for simulation
dx   = @(x,u) A*x + B*u;
dx_f = @(x,u) A*x + B*f(u);

figure(1);
clf; grid on; hold on;
bode(G)

% VRFT settings
% Reference controller (and model)
kp_ideal = 0.1;
ki_ideal = 0.1;
kd_ideal = 0.5;

M = feedback((kp_ideal + ki_ideal/s + kd_ideal*s)*G,1);

clc;
disp(evalfr(M,0))
disp(eig(M))

% Filter
W = c2d(tf(1,[0.3, 1]),Ts);

% Control strucuter
Bf = [
    tf(1,1,Ts, 'variable','z^-1');
    tf(Ts*[1 1],2*[1 -1],Ts,'variable','z^-1');
    tf([2 -2],Ts*[3 -1],Ts,'variable','z^-1')
];

%% Sample open loop data
% Additional simulation parameters
x0 = zeros(size(A,1),1);

% Input data
U = 1-2*(1-(square(2 * pi * T / 120) + 1) / 2)';

% Noise free sample data
[~,x_sample_data_linear] = ode45(@(t,x) dx(x,interp1(T,U,t)),T,x0);
y_sample_data_linear = (C*x_sample_data_linear')';

[~,x_sample_data] = ode45(@(t,x) dx_f(x,interp1(T,U,t)),T,x0);
y_sample_data = (C*x_sample_data')';

% Add noise
noise = sqrt(noise_var) * randn(size(T))';

y = y_sample_data + noise;
y_linear = y_sample_data_linear + noise;

% VRFT
[~, theta]        = VRFT1_ry(U, y,        c2d(M,Ts), Bf, W, 4, []);
[~, theta_linear] = VRFT1_ry(U, y_linear, c2d(M,Ts), Bf, W, 4, []);


disp(theta_linear./theta)
%% Figures
figure(1);
clf; grid on; hold on;
plot(T,U,'k--')
plot(T,y,'r')
plot(T,y_linear,'b')


%%
num = tf(M.num,1);
den = tf(M.den,1);

figure(2);
clf; grid on; hold on;
yline(0)
ALPHA = linspace(0,1,100);
for i = 1:length(ALPHA)
    l = zero(ALPHA(i)*den + (1-ALPHA(i))*num);

    plot(ALPHA(i),max(real(l)),'r.')
    drawnow;
end
