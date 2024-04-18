%% Descrpition
% This script is ment to demonstrate a case where performing VRFT on a
% Hammerstein system may result in instability due to underestimation of
% the real gain of the system.

%% Some general settings
clc; close;

% Text setup
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');
set(0, 'defaultAxesLabelFontSizeMultiplier', 1.3);

% Figure setup
set(0, 'defaultFigureColor', 'white');
set(0, 'defaultFigurePosition', [100, 100, 600, 400]);
set(0, 'DefaultAxesClipping', 'on');
set(0, 'DefaultLineClipping', 'on');

% Colors
NTNU_black  = [0,0,0];
NTNU_blue   = [0,80,158]/255;
NTNU_yellow = [247,208,25]/255;
NTNU_orange = [239,129,20]/255;
NTNU_brown  = [207,184,135]/255;
NTNU_purple = [176, 27, 129]/255;
NTNU_green = [188, 208, 37]/255;

%% Parameters
Ts = 0.1;
T = 0:Ts:200;
var = 0;

%% Nonlinearity
sat = @(u) min(max(u,-0.5),0.5);
f = @(u) sat(u);

%% System
% Transfer function open loop
m1 = 1;
m2 = 0.5;
c1 = 0.2;
c2 = 0.5;
k1 = 1;
k2 = 0.5;

s = tf('s');
G1 = 1/(m1*s^2+c1*s+k1);
G2 = 1/(m2*s^2+c2*s+k2);
G = minreal(G1*G2);

% State space representation
sys = ss(G);
A = sys.A;
B = sys.B;
C = sys.C;

% Controller (step response)
r = 1;
u = @(x,K) K*(r - C*x);

% Closed loop system
dx = @(x,u) A*x + B*f(u);

% Simulation of ideal controller
K_ideal = 0.09;
M = feedback((K_ideal + 0.1/s)*G,1);

clc;
disp(evalfr(M,0))
disp(eig(M))

%%
x0 = zeros(size(A,1),1);
[~,x_cl_ideal] = ode45(@(t,x) dx(x,u(x,K_ideal)),T,x0);
y_cl_ideal = (C*x_cl_ideal')';

%% Sample open loop data
% Input data
U = 2*(1-(square(2 * pi * T / 120) + 1) / 2)'-1;

% Noise
noise = sqrt(var) * randn(size(T))';

% Sample data
[~,x_sample_data] = ode45(@(t,x) dx(x,interp1(T,U,t)),T,x0);
y_sample_data = (C*x_sample_data')' + noise;

%% VRFT identification
% Filter
W = c2d(tf(1,[0.3, 1]),Ts);

% Control strucuter (P)
Bf = c2d(tf(1), Ts);

% Identification
[~, theta] = VRFT1_ry(U, y_sample_data, c2d(M,Ts), Bf, W, 4, []);

disp(theta);

% Identified controller
K_vrft = theta;

% Simulation
[~,x_cl_identified] = ode45(@(t,x) dx(x,u(x,K_vrft)),T,x0);
y_cl_identified = (C*x_cl_identified')';

%% Figures
% Plot nonlinearity
u_domain = linspace(-1,1,100);

figure(1);
clf, grid on; hold on; box on;
plot(u_domain, f(u_domain),'b')
plot(u_domain, dff(u_domain),'r')
yline(mean(df(u_domain)),'r--')
yline(1,'k-')

% Plotting sample data
figure(2);
clf; 

grid on; hold on; box on;
plot(T,y_sample_data,'b','linewidth',2,'displayname','$y(t)$')
plot(T,U,'k--','linewidth',2,'displayname','$u(t)$')

xlabel('t')
%ylabel('data')
legend('location','best')

% Plotting response
figure(3);
clf; 

subplot(2,1,1);
grid on; hold on; box on;
plot(T,y_cl_ideal,'b','linewidth',2, 'DisplayName', '$y_{\textrm{ideal}}$')
plot(T,y_cl_identified,'-','linewidth',1, 'Color', [0,0,1,0.3], 'DisplayName', '$y_{\textrm{vrft}}$')
%ylim([0,1.5])
xlabel('t')
ylabel('y')
legend('location','best')

subplot(2,1,2);
grid on; hold on; box on;
plot(T,u(x_cl_ideal',K_ideal),'r','linewidth',2, 'DisplayName', '$u_{\textrm{ideal}}$')
plot(T,u(x_cl_identified',K_vrft),'-','linewidth', 1, 'Color', [1,0,0,0.3], 'DisplayName', '$u_{\textrm{vrft}}$')
%ylim([-1,2])
xlabel('t')
ylabel('u')
legend('location','best')