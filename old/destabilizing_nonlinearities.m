addpath('VRFT_toolbox_v1\*')

clc;
clear all;

rng(1);

%% Real system
% Linear dynamics
A = [0,1;-1,-1];
B = [0,1]';
C = [1,0];

P = tf(ss(A,B,C,[]));

% Nonlinear function
f = @(u) 1./(1+exp(-u));
finv = @(v) exp(-v).^3;

% Figures
figure(3);
clf; grid on; hold on; box on;

u_range = linspace(-2,2,1000);
plot(u_range,f(u_range),'b','LineWidth',2,'DisplayName','f');
plot(u_range,finv(u_range),'r','LineWidth',2,'DisplayName','finv');
plot(u_range,f(finv(u_range)),'k--','LineWidth',2,'DisplayName','f o finv');

xlabel('u')
ylabel('y')

legend('Location','best')

%% Sample data
% Parameters
Ts = 0.1;
t = (0:Ts:500)';
var = 0.025;

% Data
u = (1-(square(2 * pi * t / 120) + 1) / 2);

% Differential equations defining the system
system_dynamics = @(time, x) A*x + B*f(interp1(t, u, time)); % interpolating u for the given time

% Simulation
x0 = [0;0];
[~,x] = ode45(system_dynamics,t,x0);

% Output
noise = sqrt(var) * randn(size(t));
y = (C*x')' + noise;

% Figure
figure(1);
clf; grid on; hold on; box on;
plot(t, u, 'k--', 'linewidth', 2)
plot(t, y, 'r', 'linewidth', 2)
xlim([0,500]);
xlabel('Time')

%% VRFT
Mr = c2d(tf(1, [3, 1]),Ts);
% W = c2d(tf(1,[0.3, 1]),Ts);
% Bf = [
%        tf(1, 1, Ts, 'variable', 'z^-1'); 
%        tf(Ts*[1 1], 2*[1 -1], Ts, 'variable', 'z^-1');
%        tf([2 -2], Ts*[3 -1], Ts, 'variable', 'z^-1');
%    ];
% 
% % VRFT
% [~, theta] = VRFT1_ry(u, y, Mr, Bf, [], 4, []);
% 
% % Extract the PID controller parameters from theta
% kp = theta(1); ki = theta(2); kd = theta(3);

%% Controller
K = [0.3,0.3,0.6];

% Simulate with controller
% Extend system dynamics
Aaug = [
    0,1,0;
    [zeros(2,1), A]];
Baug = [0; B];
Caug = [0, C];

r = [0;1;0];

% Differential equation
system_dynamics = @(time, x) Aaug*x + Baug*f(finv(K*(r-x))) - circshift(r,-1);

% Simulation
x0 = [0;0;0];
[~,x] = ode45(system_dynamics,t,x0);
y_M = (Caug*x');

% Figures
figure(2);
clf; grid on; hold on; box on;
plot(t, step(Mr,t), 'k', 'linewidth', 2)
plot(t, y_M, 'b', 'linewidth', 2)

xlim([0,50]);
xlabel('Time')