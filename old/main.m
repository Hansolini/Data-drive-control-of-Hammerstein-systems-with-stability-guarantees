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

P = 1/(s^2+s+1);

num = P.num; num = num{1};
den = P.den; den = den{1};

%% Sample data
% Parameters
Ts = 0.1;
t = (0:Ts:500)';
var = 0.025;

% Data
noise = sqrt(var) * randn(size(t));
u = (1-(square(2 * pi * t / 120) + 1) / 2);

in = [t,u];
out = sim("hammerstein_system.slx");
y = out.y + noise;

%% VRFT
Mr = c2d(tf(1, [3, 1]),Ts);
W = c2d(tf(1,[0.3, 1]),Ts);
B = [
       tf(1, 1, Ts, 'variable', 'z^-1'); 
       tf(Ts*[1 1], 2*[1 -1], Ts, 'variable', 'z^-1');
       tf([2 -2], Ts*[3 -1], Ts, 'variable', 'z^-1');
   ];

% VRFT
[~, theta] = VRFT1_ry(u, y, Mr, B, [], 4, []);

% Extract the PI controller parameters from theta
kp = theta(1); ki = theta(2); kd = theta(3);

%% Simulate with controller
with_nonlinearity = false;

t = (0:Ts:2000)';
r = ones(size(t));

%kp = 0.0683; ki = 0.1105; kd = 0.2307;
in = [t, r];
out = sim('hammerstein_system_with_PID.slx');

y_M = out.y;

%% Figures
figure(1);
clf; grid on; hold on; box on;
%plot(t, u, 'k--', 'linewidth', 2)
%plot(t, y, 'r', 'linewidth', 2)

xlim([0,500]);
xlabel('Time')

figure(2);
clf; grid on; hold on; box on;
plot(t, step(Mr,t), 'k', 'linewidth', 2)
plot(t, y_M, 'b', 'linewidth', 2)

xlim([0,50]);
xlabel('Time')


%%
x = linspace(-2,5,100);
f = @(u) u.^3 + u;
f_inv = @(v) v.^(1/3);

figure(3);
clf; grid on; hold on;
plot(x,f(f_inv(x)))