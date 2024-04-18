addpath('VRFT_toolbox_v1\*')

clc;
clear all;

rng(1);

%% Real system
m1 = 1;
m2 = 0.5;
c1 = 0.2;
c2 = 0.5;
k1 = 1;
k2 = 0.5;

s = tf('s');
P = minreal((m1*s^2+(c1+c2)*s+(k1+k2))/((m1*s^2+(c1+c2)*s+k1+k2)*(m2*s^2+c2*s+k2)-(k2+c2*s)^2));

%% Sample data
% Parameters
Ts = 0.1;
t = 0:Ts:500;
var = 0.025;

% Data
noise = sqrt(var) * randn(size(t))';
u = (1-(square(2 * pi * t / 120) + 1) / 2)';
y = lsim(P,u,t) + noise;

%% VRFT
Mr = c2d(tf(1, [3, 1]),Ts);
W = c2d(tf(1,[0.3, 1]),Ts);
B = [
        tf(1, 1, Ts, 'variable', 'z^-1'); 
        tf(Ts*[1 1], 2*[1 -1], Ts, 'variable', 'z^-1');
        tf([2 -2], Ts*[3 -1], Ts, 'variable', 'z^-1');
    ];

% Using the VRFT1_ry function, though its exact implementation is not provided in the snippet
[~, theta] = VRFT1_ry(u, y, Mr, B, W, 4, []);

% Extract the PI controller parameters from theta
Kp = theta(1); Ki = theta(2); Kd = theta(3);

% Controller
C = Kp + 1/s*Ki + s*Kd;

%% Figures
M = minreal(feedback(P*C,1));

figure(1);
clf; grid on; hold on; box on;
plot(t, u, 'k--', 'linewidth', 2)
plot(t, y, 'r', 'linewidth', 2)

xlim([0,500]);
xlabel('Time')

figure(2);
clf; grid on; hold on; box on;
plot(t, lsim(Mr,u,t), 'k', 'linewidth', 2)
plot(t, lsim(M,u,t), 'b', 'linewidth', 2)

xlim([0,240]);
xlabel('Time')